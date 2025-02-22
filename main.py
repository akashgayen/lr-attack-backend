import subprocess
import re
import os
from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel
from typing import List, Dict

app = FastAPI()

origins = [
    "http://localhost:5173",
    "http://127.0.0.1:5173",
    "https://lr-attack-frontend.vercel.app"
]

app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

class DataPayload(BaseModel):
    selectedBus: str
    busData: List[List[float]]
    lineData: List[List[float]]
    vfValues: List[float]
    PMULocation: List[int]
    suspectedNodes: List[int]

def generate_matlab_header(selectedBus: str, busData: List[List[float]], lineData: List[List[float]], vfValues: List[float], suspectedNodes: List[int], PMULocation: List[int]) -> str:
    """
    Generate MATLAB code lines (as a string) that define the matrices based on the input.
    """
    if len(busData) != 33:
        raise ValueError("busData must have 33 rows.")
    for row in busData:
        if len(row) != 3:
            raise ValueError("Each row in busData must have 3 columns.")
    if len(lineData) != 32:
        raise ValueError("lineData must have 32 rows.")
    for row in lineData:
        if len(row) != 4:
            raise ValueError("Each row in lineData must have 4 columns.")
    if len(vfValues) != 2:  
        raise ValueError("vfValues must have 33 elements.")
    if len(PMULocation) != 2:
        raise ValueError("PMULocation must have 2 elements.")
    if len(suspectedNodes) != 2:
        raise ValueError("suspectedNodes must have 2 elements.")

    matlab_lines = []
    matlab_lines.append(f"% Auto-generated header for lamda{int(selectedBus)}.m")
    matlab_lines.append("clc;")
    matlab_lines.append("clear all;")

    matlab_lines.append("M = [")
    for i, row in enumerate(busData):
        row_str = "    ".join(str(x) for x in row)
        matlab_lines.append(f"{i + 1}    {row_str};")
    matlab_lines.append("];\n")

    matlab_lines.append("l = [")
    for i, row in enumerate(lineData):
        row_str = "    ".join(str(x) for x in row)
        matlab_lines.append(f"{i + 1}    {row_str};")
    matlab_lines.append("];\n")

    matlab_lines.append(f"%%%%%%%%%%%%%%%%%%%%%Input_number3: Voltage_of_all_buses%%%%%%%%%%%%%%%%%%%%\nVOLTAGE=[1	0.997025251164255	0.982893539674410	0.975384317289895	0.967958388936322	0.949481298510629	0.945956724284158	0.932301659166856	0.925969409500574	0.920112565669725	0.919243869999172	0.917729109147354	0.911553853098794	0.909264004007971	0.907837274247153	0.906455378147553	0.904407419814540	0.903794142364604	0.996496883845413	0.992919265778768	0.992214757739055	0.991577334856689	0.979307672978918	0.972636254813258	0.969311135982705	0.947551725331661	0.944987573831078	0.933546180563927	0.925326779920730	0.921768801649497	0.917606971168094	0.916691404707977	0.916407716153128];\nvf=VOLTAGE';\nfxd1={PMULocation[0]};\nfxd2={PMULocation[1]};\nvf(fxd1)={vfValues[0]};\nvf(fxd2)={vfValues[1]};\n")
    
    matlab_lines.append(f"node1={suspectedNodes[0]};\nnode2={suspectedNodes[1]};")
    
    matlab_lines.append("% End of auto-generated header.\n")

    return "\n".join(matlab_lines)

def parse_octave_output(output: str) -> Dict[str, float]:
    """
    Parses the Octave script output and extracts relevant values into a dictionary.
    """
    parsed_data = {}

    patterns = {
        "loadability_of_node1": r"loadability_of_node1\s*=\s*([\d.]+)",
        "loadability_of_node2": r"loadability_of_node2\s*=\s*([\d.]+)",
        "actual_active_load_node1": r"Actual_activeload_of_node1\s*=\s*([\d.]+)",
        "actual_reactive_load_node1": r"Actual_reactiveload_of_node1\s*=\s*([\d.]+)",
        "actual_active_load_node2": r"Actual_activeload_of_node2\s*=\s*([\d.]+)",
        "actual_reactive_load_node2": r"Actual_reactiveload_of_node2\s*=\s*([\d.]+)",
        "lrattack_active_load_node1": r"LRattack_activeload_of_node1\s*=\s*([\d.]+)",
        "lrattack_reactive_load_node1": r"LRattack_reactiveload_of_node1\s*=\s*([\d.]+)",
        "lrattack_active_load_node2": r"LRattack_activeload_of_node2\s*=\s*([\d.]+)",
        "lrattack_reactive_load_node2": r"LRattack_reactiveload_of_node2\s*=\s*([\d.]+)",
    }

    for key, pattern in patterns.items():
        match = re.search(pattern, output)
        if match:
            parsed_data[key] = float(match.group(1))

    return parsed_data

def prepend_to_matlab_file(header_code: str, filename: str) -> None:
    """
    Prepend header_code to the specified MATLAB file.
    If the file already exists, its original content is preserved.
    """
    if os.path.exists(filename):
        with open(filename, "r") as f:
            original_content = f.read()
    else:
        original_content = ""

    with open(filename, "w") as f:
        f.write(header_code + original_content)
        
def remove_generated_header(filename: str) -> None:
    """
    Removes the auto-generated header from the MATLAB file.
    The header is identified by the markers '% Auto-generated header' and '% End of auto-generated header.'
    """
    try:
        with open(filename, 'r') as f:
            content = f.read()

        header_start = content.find("% Auto-generated header")
        header_end = content.find("% End of auto-generated header.\n")

        if header_start != -1 and header_end != -1:
            remaining_content = content[header_end + len("% End of auto-generated header.\n"):]
            
            with open(filename, 'w') as f:
                f.write(remaining_content.lstrip())
    except Exception as e:
        print(f"Error removing header from {filename}: {str(e)}")

@app.post("/get_values")
def get_values(payload: DataPayload):
    """
    This endpoint runs an Octave script and returns the extracted output.
    """
    try:
        header_code = generate_matlab_header(
            payload.selectedBus, payload.busData, payload.lineData, payload.vfValues, payload.suspectedNodes, payload.PMULocation
        )

        matlab_filename = f"lamda{payload.selectedBus}.m"

        prepend_to_matlab_file(header_code, matlab_filename)

        result = subprocess.run(["octave", "--silent", "--eval", f"run('{matlab_filename}')"],
                                capture_output=True, text=True)
        
        if result.returncode != 0:
            raise HTTPException(status_code=500, detail="Octave execution failed")

        remove_generated_header(matlab_filename)

        parsed_output = parse_octave_output(result.stdout)

        return {
            "loadability": {
                "node1": parsed_output.get("loadability_of_node1", 0),
                "node2": parsed_output.get("loadability_of_node2", 0),
            },
            "actual_load": {
                "node1": {
                    "active": parsed_output.get("actual_active_load_node1", 0),
                    "reactive": parsed_output.get("actual_reactive_load_node1", 0),
                },
                "node2": {
                    "active": parsed_output.get("actual_active_load_node2", 0),
                    "reactive": parsed_output.get("actual_reactive_load_node2", 0),
                },
            },
            "lrattack_load": {
                "node1": {
                    "active": parsed_output.get("lrattack_active_load_node1", 0),
                    "reactive": parsed_output.get("lrattack_reactive_load_node1", 0),
                },
                "node2": {
                    "active": parsed_output.get("lrattack_active_load_node2", 0),
                    "reactive": parsed_output.get("lrattack_reactive_load_node2", 0),
                },
            },
            "attack_detected": "Load distribution attack occurs between these two nodes." in result.stdout
        }

    except Exception as e:
        raise HTTPException(status_code=500, detail=str(e))


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
