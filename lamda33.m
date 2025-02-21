%%%%%%%%%%%%%%%%%%%MAIN PROGRAM%%%%%%%%%%%%%%%%%%%%%%%
lamda1=1;
lamda2=1;


no=length(M);
br=length(l);
MVAb=100;
KVb=12.66;
Zb=(KVb^2)/MVAb;

% Per unit Values
for i=1:br
    R(i,1)=(l(i,4))/Zb;
    X(i,1)=(l(i,5))/Zb;
end
for i=1:no
    P(i,1)=((M(i,2))/(1000*MVAb));
    Q(i,1)=((M(i,3))/(1000*MVAb));
end

for v=1:no
    Uj(v,1)=(abs(vf(v))).^2;
end
Uj;

BRANCH_CURRENT=[0.0392587611193980 + 0.0244295885746208i	0.0346330668854724 + 0.0222133841904172i	0.0241492770061002 + 0.0171721448718302i	0.0229166780581746 + 0.0163554634249261i	0.0222955857545815 + 0.0160480275510364i	0.0116254813764992 + 0.00560855549531235i	0.00951304033836320 + 0.00454792209167625i	0.00737256662610422 + 0.00346602719919918i	0.00672585140696803 + 0.00324638958782890i	0.00607527046462188 + 0.00302462527671843i	0.00558793584868441 + 0.00269504086911159i	0.00493664397008612 + 0.00230948207077497i	0.00428157462067728 + 0.00192024805768836i	0.00297027343792565 + 0.00102801370127445i	0.00231054530243549 + 0.000911177573515803i	0.00165102389391978 + 0.000683581944458671i	0.000990324548467306 + 0.000454562320490710i	0.00362255828177731 + 0.00161467027813899i	0.00271936886907535 + 0.00121332382217394i	0.00181339667021357 + 0.000809471930990516i	0.000906917800201573 + 0.000405027010152258i	0.00956744550334594 + 0.00463582725670506i	0.00864784883398018 + 0.00412632171823492i	0.00433054058967049 + 0.00206834208142690i	0.0100376969761370 + 0.0102303348395295i	0.00940369675118786 + 0.00996844230011973i	0.00876772003616067 + 0.00970646263200425i	0.00812386351991866 + 0.00949576295109771i	0.00682192103073245 + 0.00874819886028666i	0.00459581205130342 + 0.00225831111304047i	0.00295573197788370 + 0.00150729755796992i	0.000657602575506784 + 0.000432088790823382i];
Ibr=BRANCH_CURRENT';

%%%%%%%%%%%%%%%%%%%%%%%%% Power flow calculation %%%%%%%%%%%%%%%%%%%%%%%%
   for i=1:no-1
     P_jcf(i,1)= real((conj(vf(i+1)))*Ibr(i));
     Q_jcf(i,1)= -imag((conj(vf(i+1)))*Ibr(i));
end
  P_jcf;
  Q_jcf;
  % AAAAA
      P_kkcf = zeros(no-1,1);
      Q_kkcf = zeros(no-1,1);
   for i= 1:no-2
       P_kkcf(i,1)= P_jcf(i+1);
       Q_kkcf(i,1)=-Q_jcf(i+1);
   end
   P_kkcf;
   Q_kkcf;
   % AAAAA
   %%%%%%%%%%%%%%%%%%%%%% Power loss calculation %%%%%%%%%%%%%%%%%%%%%%%%%%

     Plosscpf = zeros(no-1,1);
     Qlosscpf = zeros(no-1,1);
     for f=1:no-2
         
         Plosscpf(f,1)=(P_kkcf(f)^2+Q_kkcf(f)^2)*R(f+1)*(Uj(f+2)^(-1));
         Qlosscpf(f,1)=(P_kkcf(f)^2+Q_kkcf(f)^2)*X(f+1)*(Uj(f+2)^(-1));
     end
Plosscpf;
Qlosscpf;
Q_jcf;
% AAAAA
 % %%%%%%%%%%%%%%%%%%%%%%%%%%% FU,FP,FQ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 for i=1:no-1
     ii=i+1;
     FU_j(i,1) =Uj(ii).^2+ ((2*(Uj(ii)))*((P_jcf(i)*R(i))+(Q_jcf(i)*X(i))-(0.5*(Uj(i)))))+((R(i).^2+X(i).^2)*(P_jcf(i).^2+Q_jcf(i).^2));
     if i==18
       FU_j(i,1) =Uj(ii).^2+ ((2*(Uj(ii)))*((P_jcf(i)*R(i))+(Q_jcf(i)*X(i))-(0.5*(Uj(2)))))+((R(i).^2+X(i).^2)*(P_jcf(i).^2+Q_jcf(i).^2));  
     end
     if i==22
       FU_j(i,1) =Uj(ii).^2+ ((2*(Uj(ii)))*((P_jcf(i)*R(i))+(Q_jcf(i)*X(i))-(0.5*(Uj(3)))))+((R(i).^2+X(i).^2)*(P_jcf(i).^2+Q_jcf(i).^2));  
     end
     if i==25
       FU_j(i,1) =Uj(ii).^2+ ((2*(Uj(ii)))*((P_jcf(i)*R(i))+(Q_jcf(i)*X(i))-(0.5*(Uj(6)))))+((R(i).^2+X(i).^2)*(P_jcf(i).^2+Q_jcf(i).^2));  
     end
 end
 Uj;
 FU_j;
 % AAAA
 for i=1:no-1

     if i==node1-1
         FP_j(i,1)= P_jcf(i)-Plosscpf(i)-P_kkcf(i)-(lamda1*P(node1));
     elseif i==node2-1
         FP_j(i,1)= P_jcf(i)-Plosscpf(i)-P_kkcf(i)-(lamda2*P(node2));
     % elseif i==node3
     %     FP_j(i,1)= P_jcf(i)-Plosscpf(i)-P_kkcf(i)-(lamda3*P(node3));
     elseif i==1
         FP_j(i,1)= P_jcf(i)-Plosscpf(i)-Plosscpf(17)-P_jcf(i+1)-P_jcf(18);
     elseif i==2
         FP_j(i,1)= P_jcf(i)-Plosscpf(i)-Plosscpf(21)-P_jcf(i+1)-P_jcf(22);
     elseif i==5
         FP_j(i,1)= P_jcf(i)-Plosscpf(i)-Plosscpf(24)-P_jcf(i+1)-P_jcf(25);
     elseif i==17
         FP_j(i,1)= P_jcf(i);
     elseif i==24
         FP_j(i,1)= P_jcf(i);
     elseif i==21
         FP_j(i,1)= P_jcf(i);

     else
         FP_j(i,1)= P_jcf(i)-Plosscpf(i)-P_kkcf(i);
     end
 end
 AOA=[FP_j,P(2:no)];
 % AAAAA
 for i=1:no-1

     if i==node1-1
         FQ_j(i,1)= Q_jcf(i)-Qlosscpf(i)+Q_kkcf(i)-(lamda1*Q(node1));
     elseif i==node2-1
         FQ_j(i,1)= Q_jcf(i)-Qlosscpf(i)+Q_kkcf(i)-(lamda2*Q(node2));
     % elseif i==node3
     %     FQ_j(i,1)= Q_jcf(i)-Qlosscpf(i)+Q_kkcf(i)-(lamda3*Q(node3));
     elseif i==1
         FQ_j(i,1)= Q_jcf(i)-Qlosscpf(i)-Qlosscpf(17)-Q_jcf(i+1)-Q_jcf(18);
     elseif i==2
         FQ_j(i,1)= Q_jcf(i)-Qlosscpf(i)-Qlosscpf(21)-Q_jcf(i+1)-Q_jcf(22);
     elseif i==5
         FQ_j(i,1)= Q_jcf(i)-Qlosscpf(i)-Qlosscpf(24)-Q_jcf(i+1)-Q_jcf(25);
     elseif i==17
         FQ_j(i,1)= Q_jcf(i);
     elseif i==24
         FQ_j(i,1)= Q_jcf(i);
     elseif i==21
         FQ_j(i,1)= Q_jcf(i);

     else
         FQ_j(i,1)= Q_jcf(i)-Qlosscpf(i)+Q_kkcf(i);
     end
     end
 
 FQ_j;

 COA=[FQ_j,Q(2:no)];
P_jcf;
Q_jcf;
% AAAAAA
itertion=1;
itermaximum=20000;
epsilon2=1e-6;
maxraz2=epsilon2+1;
% vbnewcf2=vf;
while (itertion<itermaximum)&&(maxraz2>epsilon2)

FS=[FU_j;FP_j;FQ_j];
JO=zeros((no-1),1);
POWW=[JO;P(2:no);Q(2:no)];
missmatch=(POWW-FS);
% AAAA

  %%%%%%%%%%%%%%%%%%%%%%% Corrector Step %%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%%%% Jacobian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  %%%%%%%%%%%%%%%%%%%%%%% delFU/delU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  J1=zeros(no-1,no-2);
  for i=1:(no-1)
      for k=1:(no-2)
          if (i==18) && (k==17)
             J1(i,k)=(2*Uj(k+2))+(2*((P_jcf(k+1)*R(k+1))+(Q_jcf(k+1)*X(k+1))-(0.5*Uj(2))));
          end
          if (i==18) && (k==1)
            J1(i,k)=-Uj(19);
          end
          if (i==25) && (k==5)
             J1(i,k)=-Uj(26);
          end
          if (i==22) && (k==2)
            J1(i,k)=-Uj(23);
         end
          if k==i && i<17 && k<17
              J1(i,k)=(2*Uj(k+1))+(2*((P_jcf(k)*R(k))+(Q_jcf(k)*X(k))-(0.5*Uj(k))));
          elseif k==i-1 && i<18 && k<17
              J1(i,k)=-Uj(k+2);
          elseif k==i-1 && i>18 && k>17
              J1(i,k)=(2*Uj(k+2))+(2*((P_jcf(k+1)*R(k+1))+(Q_jcf(k+1)*X(k+1))-(0.5*Uj(k+1))));
              if (i==22) && (k==21)
                  J1(i,k)=(2*Uj(k+2))+(2*((P_jcf(k+1)*R(k+1))+(Q_jcf(k+1)*X(k+1))-(0.5*Uj(3))));
              end

              if (i==25) && (k==24)
                  J1(i,k)=(2*Uj(k+2))+(2*((P_jcf(k+1)*R(k+1))+(Q_jcf(k+1)*X(k+1))-(0.5*Uj(6))));
              end
          elseif k==i-2 && i>18 && k>16
              J1(i,k)=-Uj(k+3);
               if (i==25) && (k==23)
                 J1(i,k)=0;
               end
               if (i==22) && (k==20)
                 J1(i,k)=0;
              end 
         end
      end
  end
  
  J1;
  P_jcf;
  % AAAAAAA

   %%%%%%%%%%%%%%%%%%%%% delFU/delP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   J2=zeros(no-1,no-1);
    for i=1:(no-1)
       for k=1:(no-1)
           if k==i
        J2(i,k)=J2(i,k)+(2*R(k)*Uj(k+1))+(2*P_jcf(k)*(R(k).^2+X(k).^2));
          end
       end
    end
    J2;
    % AAAA

   %%%%%%%%%%%%%%%%%%%%% delFU/delQ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   J3=zeros(no-1,no-1);
    for i=1:(no-1)
       for k=1:(no-1)
           if k==i
           J3(i,k)=J3(i,k)+(2*X(k)*Uj(k+1))+(2*Q_jcf(k)*(R(k).^2+X(k)^2));
           end
       end
    end
   J3;
   % AAAA
    %%%%%%%%%%%%%%%%%%%%% delFP/delU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    J44=zeros(no-1,no-2);
    for i=1:(no-1)
        for k=1:(no-2)

            if i==1 && k==17
                J44(i,k)=R(k+1)*(P_kkcf(k).^2+Q_kkcf(k).^2)*Uj(k+2).^(-2);
            end
            if i==2 && k==21
                J44(i,k)=R(k+1)*(P_kkcf(k).^2+Q_kkcf(k).^2)*Uj(k+2).^(-2);
            end
            if i==5 && k==24
                J44(i,k)=R(k+1)*(P_kkcf(k).^2+Q_kkcf(k).^2)*Uj(k+2).^(-2);
            end
            if k==i && i<17 && k<17
                % J4(i,k)=J4(i,k)-Gm(k);
                J44(i,k)=J44(i,k)-0;
            elseif k==i && i>16 && k>16
                J44(i,k)=R(k+1)*(P_kkcf(k).^2+Q_kkcf(k).^2)*Uj(k+2).^(-2);
            elseif k==(i+1) && i<16 && k<17
                J44(i,k)=R(k)*(P_kkcf(k-1).^2+Q_kkcf(k-1).^2)*Uj(k+1).^(-2);

            end

        end
    end

% Gm;
J44;
J4=J44;

 %%%%%%%%%%%%%%%%%%%%% delFP/delP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   J55=zeros(no-1,no-1);
   for i=1:(no-1)
       for k=1:(no-1)
           if i==1 && k==18
                  J55(i,k)=-1-(2*R(k)*P_kkcf(k-1)*Uj(k+1).^(-1));
               end
               if i==2 && k==22
                  J55(i,k)=-1-(2*R(k)*P_kkcf(k-1)*Uj(k+1).^(-1));
               end
               if i==5 && k==25
                  J55(i,k)=-1-(2*R(k)*P_kkcf(k-1)*Uj(k+1).^(-1));
               end
           if k==i
           J55(i,k)=J55(i,k)+1;
           elseif k==i+1
               J55(i,k)=-1-(2*R(k)*P_kkcf(k-1)*Uj(k+1).^(-1));
               
           
           end
       end
   end
   J55;
   J5=J55;
 %%%%%%%%%%%%%%%%%%%%% delFP/delQ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   J66=zeros(no-1,no-1);
   for i=1:(no-1)
       for k=1:(no-1)
           if i==1 && k==18
                  J66(i,k)=(2*R(k)*Q_kkcf(k-1)*Uj(k+1).^(-1));
                  end
                  if i==2 && k==22
                  J66(i,k)=(2*R(k)*Q_kkcf(k-1)*Uj(k+1).^(-1));
                  end
                  if i==5 && k==25
                  J66(i,k)=(2*R(k)*Q_kkcf(k-1)*Uj(k+1).^(-1));
                  end
           if k==i
           J66(i,k)=J66(i,k)+0;
           elseif k==i+1
               % J6(i,k)=-(2*R(k)*Q_kkcf(k-1)*Uj(k+1).^(-1));
               J66(i,k)=(2*R(k)*Q_kkcf(k-1)*Uj(k+1).^(-1));
                  
           end
       end
   end
   J66;
   J6=J66;
   %%%%%%%%%%%%%%%%%%%%% delFQ/delU %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   J77=zeros(no-1,no-2);
   for i=1:(no-1)
       for k=1:(no-2)
           if i==1 && k==17
               J77(i,k)=X(k+1)*(P_kkcf(k).^2+Q_kkcf(k).^2)*Uj(k+2).^(-2);
           end
           if i==2 && k==21
               J77(i,k)=X(k+1)*(P_kkcf(k).^2+Q_kkcf(k).^2)*Uj(k+2).^(-2);
           end
           if i==5 && k==24
               J77(i,k)=X(k+1)*(P_kkcf(k).^2+Q_kkcf(k).^2)*Uj(k+2).^(-2);
           end

           if k==i && i<17 && k<17
               % J4(i,k)=J4(i,k)-Gm(k);
               J77(i,k)=J77(i,k)-0;
           elseif k==i && i>16 && k>16
               J77(i,k)=X(k+1)*(P_kkcf(k).^2+Q_kkcf(k).^2)*Uj(k+2).^(-2);
           elseif k==(i+1) && i<16 && k<17
               J77(i,k)=X(k)*(P_kkcf(k-1).^2+Q_kkcf(k-1).^2)*Uj(k+1).^(-2);

           end
       end
   end
   J77;
   J7=J77;
   %%%%%%%%%%%%%%%%%%%%% delFQ/delP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   J88=zeros(no-1,no-1);
   for i=1:(no-1)
       for k=1:(no-1)
           if i==1 && k==18 
                      J88(i,k)=-(2*X(k)*P_kkcf(k-1)*Uj(k+1).^(-1));
                    end
                    if i==2 && k==22 
                      J88(i,k)=-(2*X(k)*P_kkcf(k-1)*Uj(k+1).^(-1));
                    end
                    if i==5 && k==25 
                      J88(i,k)=-(2*X(k)*P_kkcf(k-1)*Uj(k+1).^(-1));
                    end
           if k==i
           J88(i,k)=J88(i,k)+0;
           elseif k==i+1
               J88(i,k)=-(2*X(k)*P_kkcf(k-1)*Uj(k+1).^(-1));
               
           end
       end
   end
   J88;
  
     J8=J88;
  
   %%%%%%%%%%%%%%%%%%%%% delFQ/delQ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   J99=zeros(no-1,no-1);
   for i=1:(no-1)
       for k=1:(no-1)
           if i==1 && k==18
               J99(i,k)=-1+(2*X(k)*Q_kkcf(k-1)*Uj(k+1).^(-1));
           end
           if i==2 && k==22
               J99(i,k)=-1+(2*X(k)*Q_kkcf(k-1)*Uj(k+1).^(-1));
           end
           if i==5 && k==25
               J99(i,k)=-1+(2*X(k)*Q_kkcf(k-1)*Uj(k+1).^(-1));
           end
           if k==i
               J99(i,k)=J99(i,k)+1;
           elseif k==i+1
               % J9(i,k)=-1-(2*X(k)*Q_kkcf(k-1)*Uj(k+1).^(-1));
               J99(i,k)=-1+(2*X(k)*Q_kkcf(k-1)*Uj(k+1).^(-1));

           end
       end
   end
   J99;
  
     J9=J99;
  
    JB=[J1 J2 J3; J4 J5 J6; J7 J8 J9];
    JB(:,31)=[];
    size(JB);
% AAAAA
Umat=zeros(no-1,2);
FPmat=zeros(no-1,2);
 for i=1:no-1
     for h=1:2
         if i==node1-1 && h==1
           FPmat(i,h)=-P(node1); 
         elseif i==node2-1 && h==2
             FPmat(i,h)=-P(node2);
         %  elseif i==node3 && h==3
         %     FPmat(i,h)=-P(node3);
      end 
     end
 end
     FPmat;

     FQmat=zeros(no-1,2);
 for i=1:no-1
     for h=1:2
         if i==node1-1 && h==1
           FQmat(i,h)=-Q(node1); 
         elseif i==node2-1 && h==2
             FQmat(i,h)=-Q(node2);
         %  elseif i==node3 && h==3
         %  FQmat(i,h)=-Q(node3);
      end 
     end
 end
     FQmat;

MEN=[Umat;FPmat;FQmat];

   % AAAAA
       
       Jacobian_matrix=[JB,MEN];
       size(Jacobian_matrix);
       % AAAA
       inverse_Jacobian_matrix=inv(Jacobian_matrix);
       solution=inverse_Jacobian_matrix*missmatch;
       % AAAAAA

 %%%%%%%%%%%%%%%%%%%%% Update part %%%%%%%%%%%%%%%%%%%%%%%
Uj;
SP=[Uj(2:17);Uj(19:22);Uj(23:32);P_jcf;Q_jcf;lamda1;lamda2];
size(SP);
% AAAAA
etaa=0.1;
final_matrix=SP+(etaa.*solution);
% Uj_nw=final_matrix(1:15,1);
Uj_new=[1;final_matrix(1:16,1);(0.9)^2;final_matrix(17:19,1);final_matrix(20:30,1);(0.91)^2];
% AAAAA
P_j_new=final_matrix(31:62,1);
Q_j_new=final_matrix(63:94,1);
lamdanew1=final_matrix(95,1);
lamdanew2=final_matrix(96,1);
% lamdanew3=final_matrix(98,1);
lamda1=lamdanew1;
lamda2=lamdanew2;
% lamda3=lamdanew3;
vbnewcf=[sqrt(Uj_new)];

Uj=Uj_new;
P_jcf=P_j_new;
Q_jcf=Q_j_new;


for i=2:no
    if i==2
        delta(i,1)=0-(asind(((P_jcf(i-1)*X(i-1))-(Q_jcf(i-1)*R(i-1)))/((vbnewcf(i))*(vbnewcf(i-1)))));
    elseif i==19
        delta(i,1)= delta(2,1)-(asind(((P_jcf(i-1)*X(i-1))-(Q_jcf(i-1)*R(i-1)))/((vbnewcf(i))*(vbnewcf(i-1)))));
        elseif i==23
        delta(i,1)= delta(3,1)-(asind(((P_jcf(i-1)*X(i-1))-(Q_jcf(i-1)*R(i-1)))/((vbnewcf(i))*(vbnewcf(i-1)))));
        elseif i==26
        delta(i,1)= delta(6,1)-(asind(((P_jcf(i-1)*X(i-1))-(Q_jcf(i-1)*R(i-1)))/((vbnewcf(i))*(vbnewcf(i-1)))));
    else
        delta(i,1)= delta(i-1,1)-(asind(((P_jcf(i-1)*X(i-1))-(Q_jcf(i-1)*R(i-1)))/((vbnewcf(i))*(vbnewcf(i-1)))));
    end
end
delta;

% AAAAA
% AAAA
for i=1:no
    vbnewcf2(i,1)=(vbnewcf(i)*cosd(delta(i)))+(1j*vbnewcf(i)*sind(delta(i)));
end

vf=vbnewcf2;

%%%%%%%%%%%%%%%%%%%%%%%%% Power flow calculation %%%%%%%%%%%%%%%%%%%%%%%%

     P_kkcfnew = zeros(no-1,1);
     Q_kkcfnew = zeros(no-1,1);
  
   for i= 1:no-2
       P_kkcfnew(i,1)= P_jcf(i+1);
       Q_kkcfnew(i,1)=-Q_jcf(i+1);
   end
   P_kkcf=P_kkcfnew;
   Q_kkcf=Q_kkcfnew;
   %%%%%%%%%%%%%%%%%%%%%% Power loss calculation %%%%%%%%%%%%%%%%%%%%%%%%%%

   Plosscpfnew = zeros(no-1,1);
   Qlosscpfnew = zeros(no-1,1);
   for f=1:no-2
       Plosscpfnew(f,1)=(P_kkcf(f)^2+Q_kkcf(f)^2)*R(f+1)*(Uj(f+2)^(-1));
       Qlosscpfnew(f,1)=(P_kkcf(f)^2+Q_kkcf(f)^2)*X(f+1)*(Uj(f+2)^(-1));
   end
   Plosscpf=Plosscpfnew;
   Qlosscpf=Qlosscpfnew;

 % %%%%%%%%%%%%%%%%%%%%%%%%%%% FU,FP,FQ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 % %%%%%%%%%%%%%%%%%%%%%%%%%%% FU,FP,FQ %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:no-1
     ii=i+1;
     FU_jnew(i,1) =Uj(ii).^2+ ((2*(Uj(ii)))*((P_jcf(i)*R(i))+(Q_jcf(i)*X(i))-(0.5*(Uj(i)))))+((R(i).^2+X(i).^2)*(P_jcf(i).^2+Q_jcf(i).^2));
     if i==18
       FU_jnew(i,1) =Uj(ii).^2+ ((2*(Uj(ii)))*((P_jcf(i)*R(i))+(Q_jcf(i)*X(i))-(0.5*(Uj(2)))))+((R(i).^2+X(i).^2)*(P_jcf(i).^2+Q_jcf(i).^2));  
     end
     if i==22
       FU_jnew(i,1) =Uj(ii).^2+ ((2*(Uj(ii)))*((P_jcf(i)*R(i))+(Q_jcf(i)*X(i))-(0.5*(Uj(3)))))+((R(i).^2+X(i).^2)*(P_jcf(i).^2+Q_jcf(i).^2));  
     end
     if i==25
       FU_jnew(i,1) =Uj(ii).^2+ ((2*(Uj(ii)))*((P_jcf(i)*R(i))+(Q_jcf(i)*X(i))-(0.5*(Uj(6)))))+((R(i).^2+X(i).^2)*(P_jcf(i).^2+Q_jcf(i).^2));  
     end
 end
 Uj;
 FU_jnew;
 % AAAA
 for i=1:no-1

     if i==node1-1
         FP_jnew(i,1)= P_jcf(i)-Plosscpf(i)-P_kkcf(i)-(lamda1*P(node1));
     elseif i==node2-1
         FP_jnew(i,1)= P_jcf(i)-Plosscpf(i)-P_kkcf(i)-(lamda2*P(node2));
     % elseif i==node3
     %     FP_jnew(i,1)= P_jcf(i)-Plosscpf(i)-P_kkcf(i)-(lamda3*P(node3));
     elseif i==1
         FP_jnew(i,1)= P_jcf(i)-Plosscpf(i)-Plosscpf(17)-P_jcf(i+1)-P_jcf(18);
     elseif i==2
         FP_jnew(i,1)= P_jcf(i)-Plosscpf(i)-Plosscpf(21)-P_jcf(i+1)-P_jcf(22);
     elseif i==5
         FP_jnew(i,1)= P_jcf(i)-Plosscpf(i)-Plosscpf(24)-P_jcf(i+1)-P_jcf(25);
     elseif i==17
         FP_jnew(i,1)= P_jcf(i);
     elseif i==24
         FP_jnew(i,1)= P_jcf(i);
     elseif i==21
         FP_jnew(i,1)= P_jcf(i);

     else
         FP_jnew(i,1)= P_jcf(i)-Plosscpf(i)-P_kkcf(i);
     end
 end
 AOA=[FP_jnew,P(2:no)];
 % AAAAA
 for i=1:no-1

     if i==node1-1
         FQ_jnew(i,1)= Q_jcf(i)-Qlosscpf(i)+Q_kkcf(i)-(lamda1*Q(node1));
     elseif i==node2-1
         FQ_jnew(i,1)= Q_jcf(i)-Qlosscpf(i)+Q_kkcf(i)-(lamda2*Q(node2));
     % elseif i==node3
     %     FQ_jnew(i,1)= Q_jcf(i)-Qlosscpf(i)+Q_kkcf(i)-(lamda3*Q(node3));
     elseif i==1
         FQ_jnew(i,1)= Q_jcf(i)-Qlosscpf(i)-Qlosscpf(17)-Q_jcf(i+1)-Q_jcf(18);
     elseif i==2
         FQ_jnew(i,1)= Q_jcf(i)-Qlosscpf(i)-Qlosscpf(21)-Q_jcf(i+1)-Q_jcf(22);
     elseif i==5
         FQ_jnew(i,1)= Q_jcf(i)-Qlosscpf(i)-Qlosscpf(24)-Q_jcf(i+1)-Q_jcf(25);
     elseif i==17
         FQ_jnew(i,1)= Q_jcf(i);
     elseif i==24
         FQ_jnew(i,1)= Q_jcf(i);
     elseif i==21
         FQ_jnew(i,1)= Q_jcf(i);

     else
         FQ_jnew(i,1)= Q_jcf(i)-Qlosscpf(i)+Q_kkcf(i);
     end
 end
 FQ_jnew;

   FU_j=FU_jnew;
   FP_j=FP_jnew;
   FQ_j=FQ_jnew;
% AAAA
   V_check=max(abs(vbnewcf2-vf));
   
   P_check=max(abs(P(2:no)-FP_jnew));
   
   Q_check=max(abs(Q(2:no)-FQ_jnew));
  
   maxraz2=max([V_check P_check Q_check]);


   itertion=itertion+1;
   
end


%%%%%%%%%%%%%Output no1: loadability value%%%%%%%%%%%%%%%%%%%%%%%
loadability_of_node1=lamda1
loadability_of_node2=lamda2

% Check conditions and display the result
if (lamda1 > 1 && lamda1 > 1) || (lamda1 < 1 && lamda1 < 1)
    disp('Load distribution attack occurs between these two nodes.');
elseif lamda1 == 1 && lamda1 == 1
    disp('No Load distribution attack occurs between these two nodes.');
else
    disp('Condition does not indicate a load distribution attack.');
end


%%%%%%%%%%%%%%%%%% Output no3: Actual loading and Load Redistribution Attack loading%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Actual_activeload_of_node1=(P(node1)*MVAb)/(1+lamda1)
Actual_reactiveload_of_node1=(Q(node1)*MVAb)/(1+lamda1)

Actual_activeload_of_node2=(P(node2)*MVAb)/(1+lamda2)
Actual_reactiveload_of_node2=(Q(node2)*MVAb)/(1+lamda2)

LRattack_activeload_of_node1=P(node1)*MVAb
LRattack_reactiveload_of_node1=Q(node1)*MVAb

LRattack_activeload_of_node2=P(node2)*MVAb
LRattack_reactiveload_of_node2=Q(node2)*MVAb