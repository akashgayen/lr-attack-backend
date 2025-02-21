FROM python:3.9-slim

ENV PORT=8000
WORKDIR /app

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    octave \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY . .
RUN chmod +x *.m

RUN echo '#!/bin/bash\n\
    uvicorn main:app --host 0.0.0.0 --port $PORT' > start.sh && \
    chmod +x start.sh

CMD ["./start.sh"]