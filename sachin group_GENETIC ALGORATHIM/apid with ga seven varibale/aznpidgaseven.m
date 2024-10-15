%%% genetic algotihm based adaptive PID controller 
%%%%plant model is exp(-0.2s) /(s+1)^2  
%%%%%%GA-APID [here we are calculating all seven variables with in the defined limit using genetic algorithm kp, ki, kd  k1, k2, k3 and k4  
%%%GA-APID


clc ;

clear;

popsize=10;

Nt=28

no_of_variable=7

pop=round(rand(popsize,Nt));

aa=pop(:,1:4);

bb=pop(:,5:8); 

cc=pop(:,9:12);




dd=pop(:,13:16);

ee=pop(:,17:20); 

ff=pop(:,21:24);

gg=pop(:,25:28);




d1=bi2de(aa,'left-msb');
d2=bi2de(bb,'left-msb');
d3=bi2de(cc,'left-msb');


d4=bi2de(dd,'left-msb');
d5=bi2de(ee,'left-msb');
d6=bi2de(ff,'left-msb');
d7=bi2de(gg,'left-msb');

                  
%%%%%%%%%%%%%%%%%%%%%initial value


 ku=10.5;
 tu=2.0333;
 

kp=0.6*ku


 ti=tu/2;
 td=tu/8;
 
ki=kp*(0.1/ti);


kd=kp*(td/0.1);
    
                  
                  
                

xh1=kp+0.20*kp;
xl1=kp-0.20*kp;


xh2=ki+0.20*ki;
xl2=ki-0.20*ki;

xh3=kd+0.20*kd;
xl3=kd-0.20*kd;





xh4=5;
xl4=0;

xh5=5;
xl5=0


xh6=30;
xl6=0;



xh7=1;
xl7=0.1;


nmut=3

delay=0.2;


 h=0.1; 
 t = 0:h:16;  
 tf=160;
 
%%%%%%%%%%%%%%%%%%5

for i=1:10
       

x1(i,1)=xl1+((xh1-xl1)/(2^4-1))*d1(i);

x2(i,1)=xl2+((xh2-xl2)/(2^4-1))*d2(i);

x3(i,1)=xl3+((xh3-xl3)/(2^4-1))*d3(i);





x4(i,1)=xl4+((xh4-xl4)/(2^4-1))*d4(i);

x5(i,1)=xl5+((xh5-xl5)/(2^4-1))*d5(i);

x6(i,1)=xl6+((xh6-xl6)/(2^4-1))*d6(i);

x7(i,1)=xl7+((xh7-xl7)/(2^4-1))*d7(i);



end




d=[x1 x2 x3 x4 x5 x6 x7] ;




[iae] = objecsevenvariable(d)




[iae ind]=sort(iae);

pop=pop(ind,:);


d=d(ind,:)


iga=0;
 
 
 while iga<60
    
    
    
iga=iga+1;

%%%%%%%%%%%%%%%%selection of mimimum five iae%%%%%%%%%%%%%%%%%%%%%%

  
  pop;

  p1=[pop(1,1:2) pop(2,3:4) pop(1,5:6) pop(2,7:8) pop(1,9:10) pop(2,11:12)   pop(1,13:14)   pop(2,15:16)  pop(1,17:18)  pop(2,19:20)  pop(1,21:22) pop(2,23:24)  pop(1,25:26)  pop(2,27:28) ];        
 
  p2=[pop(2,1:2) pop(1,3:4) pop(2,5:6) pop(1,7:8) pop(2,9:10) pop(1,11:12)   pop(2,13:14)   pop(1,15:16)  pop(2,17:18)  pop(1,19:20)  pop(2,21:22) pop(1,23:24)  pop(2,25:26)  pop(1,27:28) ];   
 
  p3=[pop(3,1:2) pop(4,3:4) pop(3,5:6) pop(4,7:8) pop(3,9:10) pop(4,11:12)    pop(3,13:14)   pop(4,15:16)  pop(3,17:18)  pop(4,19:20)  pop(3,21:22) pop(4,23:24)  pop(3,25:26)  pop(4,27:28) ];        
  
  p4=[pop(4,1:2) pop(3,3:4) pop(4,5:6) pop(3,7:8) pop(4,9:10) pop(3,11:12)   pop(4,13:14)   pop(3,15:16)  pop(4,17:18)  pop(3,19:20)  pop(4,21:22) pop(3,23:24)  pop(4,25:26)  pop(3,27:28) ];        
  
  p5=[pop(5,1:2) pop(1,3:4) pop(5,5:6) pop(1,7:8) pop(5,9:10) pop(1,11:12)    pop(5,13:14)   pop(1,15:16)  pop(5,17:18)  pop(1,19:20)  pop(5,21:22) pop(1,23:24)  pop(5,25:26)  pop(1,27:28) ];        
  
  
  
  
  
  m =[p1;p2;p3;p4;p5];
  
  pop(6:10,:)=m ;
  



aa=pop(:,1:4);

bb=pop(:,5:8); 

cc=pop(:,9:12);




dd=pop(:,13:16);

ee=pop(:,17:20); 

ff=pop(:,21:24);

gg=pop(:,25:28);




d1=bi2de(aa,'left-msb');
d2=bi2de(bb,'left-msb');
d3=bi2de(cc,'left-msb');


d4=bi2de(dd,'left-msb');
d5=bi2de(ee,'left-msb');
d6=bi2de(ff,'left-msb');
d7=bi2de(gg,'left-msb');






for i=1:10
       

x1(i,1)=xl1+((xh1-xl1)/(2^4-1))*d1(i);

x2(i,1)=xl2+((xh2-xl2)/(2^4-1))*d2(i);

x3(i,1)=xl3+((xh3-xl3)/(2^4-1))*d3(i);


x4(i,1)=xl4+((xh4-xl4)/(2^4-1))*d4(i);

x5(i,1)=xl5+((xh5-xl5)/(2^4-1))*d5(i);

x6(i,1)=xl6+((xh6-xl6)/(2^4-1))*d6(i);

x7(i,1)=xl7+((xh7-xl7)/(2^4-1))*d7(i);



end




d=[x1 x2 x3 x4 x5 x6 x7] ;



 %%%%objective function
 
 
[iae] = objecsevenvariable(d)


[iae ind]=sort(iae);

pop=pop(ind,:);

d=d(ind,:)




%%%%%%%%%%%%%%%%%%%%%%%%%%%mutation%%%%%%%%%%%%%%


pop ;


mrow=ceil(rand(1,nmut)*(popsize-1))+1 ;

mcol=ceil(rand(1,nmut)*Nt) ;


for ii=1:nmut
    
pop(mrow(ii),mcol(ii))=abs(pop(mrow(ii),mcol(ii))-1);

end


pop;


aa=pop(:,1:4);

bb=pop(:,5:8); 

cc=pop(:,9:12);


dd=pop(:,13:16);

ee=pop(:,17:20); 

ff=pop(:,21:24);

gg=pop(:,25:28);




d1=bi2de(aa,'left-msb');
d2=bi2de(bb,'left-msb');
d3=bi2de(cc,'left-msb');


d4=bi2de(dd,'left-msb');
d5=bi2de(ee,'left-msb');
d6=bi2de(ff,'left-msb');
d7=bi2de(gg,'left-msb');






for i=1:10
       

x1(i,1)=xl1+((xh1-xl1)/(2^4-1))*d1(i);

x2(i,1)=xl2+((xh2-xl2)/(2^4-1))*d2(i);

x3(i,1)=xl3+((xh3-xl3)/(2^4-1))*d3(i);


x4(i,1)=xl4+((xh4-xl4)/(2^4-1))*d4(i);

x5(i,1)=xl5+((xh5-xl5)/(2^4-1))*d5(i);

x6(i,1)=xl6+((xh6-xl6)/(2^4-1))*d6(i);

x7(i,1)=xl7+((xh7-xl7)/(2^4-1))*d7(i);


end




d=[x1 x2 x3 x4 x5 x6 x7] ;



 %%%%objective function
 

[iae] = objecsevenvariable(d)


[iae ind]=sort(iae);

pop=pop(ind,:);

d=d(ind,:)




end 


ind;

iga

d;



% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  for j=1  
    
    
    kp=d(j,1);
    ki=d(j,2);
    kd=d(j,3);
    
    
    
    k1 =d(j,4);
    k2=d(j,5);
    k3=d(j,6);
    k4=d(j,7)
    


 y = zeros(1,length(t)); 
 u = zeros(1,length(t));
 u = zeros(1,length(t));
 
 kpp = zeros(1,length(t)); 
 kii = zeros(1,length(t));
 kdd = zeros(1,length(t));
 
 v = zeros(1,length(t));
 x = zeros(1,length(t));

r= 1;

setpoint=r;
y(1) = 0;
x(1)=0 ;
e(1)=r - y(1);


u(1)=kp*e(1)+ki*sum(e);

 kpp(1)=0;
 kii(1)=0;
 kdd(1)=0;
 v(1)=0;


F_xy = @(x) -2*x; 

for i = 1:2
         
    
        k_1 = F_xy(x(i));
        k_2 = F_xy(x(i)+0.5*h*k_1);
        k_3 = F_xy((x(i)+0.5*h*k_2));
        k_4 = F_xy((x(i)+k_3*h));
        
        x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
  
        y(i+1)=y(i)+h*x(i+1);
    
        e(i+1)=r-y(i+1);
   
        er=sum(e);
   
        ed=e(i+1)-e(i);
     

        v(i+1)=(e(i+1)/r)*(ed/r);
        
      

kpp(i+1)=kp*(1+k1*abs(v(i+1)));
kii(i+1)=ki*(k4+k2*v(i+1));
kdd(i+1)=kd*(1+k3*abs(v(i+1)));


u(i+1)=kpp(i+1)*e(i+1)+ kii(i+1)*sum(e)+ kdd(i+1)*ed;

end




F_xy = @(u,y,x) u-y-2*x;


for i = 3:(tf/2)
    
    
    
    k_1 = F_xy(u(i-2),y(i),x(i));
    
    k_2 = F_xy(u(i-2)+0.5*h,y(i)+0.5*h,x(i)+0.5*h*k_1);
    
    k_3 = F_xy((u(i-2)+0.5*h),(y(i)+0.5*h),(x(i)+0.5*h*k_2));
    
    k_4 = F_xy((u(i-2)+h ) ,(y(i)+h) ,(x(i)+k_3*h));
    
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
   y(i+1)=y(i)+h*x(i+1); 

  
e(i+1)=r-y(i+1);
ed=e(i+1)-e(i);  
er=sum(e);
v(i+1)=(e(i+1)/r)*(ed/r);
kpp(i+1)=kp*(1+k1*abs(v(i+1)));

kii(i+1)=ki*(k4+k2*v(i+1));

kdd(i+1)=kd*(1+k3*abs(v(i+1)));


u(i+1)=kpp(i+1)*e(i+1)+ kii(i+1)*sum(e)+ kdd(i+1)*ed;
    
end 


u((tf/2)+1)= -20;

F_xy = @(u,y,x) u-y-2*x;


for i = (tf/2)+1:length(t)
    
    
    k_1 = F_xy(u(i-2),y(i),x(i));
    
    k_2 = F_xy(u(i-2)+0.5*h,y(i)+0.5*h,x(i)+0.5*h*k_1);
    
    k_3 = F_xy((u(i-2)+0.5*h),(y(i)+0.5*h),(x(i)+0.5*h*k_2));
    
    k_4 = F_xy((u(i-2)+h ) ,(y(i)+h) ,(x(i)+k_3*h));
    
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
   y(i+1)=y(i)+h*x(i+1); 

  
e(i+1)=r-y(i+1);
ed=e(i+1)-e(i);  
er=sum(e);
v(i+1)=(e(i+1)/r)*(ed/r);
kpp(i+1)=kp*(1+k1*abs(v(i+1)));

kii(i+1)=ki*(k4+k2*v(i+1));

kdd(i+1)=kd*(1+k3*abs(v(i+1)));


u(i+1)=kpp(i+1)*e(i+1)+ kii(i+1)*sum(e)+ kdd(i+1)*ed;
    
end 


z=y(1:(tf+1));

figure(1)

plot(t,z) 

xlabel('Time t ')
ylabel('Response y')
title(' exp(-0.2s) /(s+1)^2  GA-APID [GA based adaptive PID, variation in seven variables i.e.kp, ki, kd  k1, k2, k3 and k4  ')
grid on

  end

 
  
%%%%%iae/itae

y;

yy = y(1:(tf/2));

[yymax tp]=max(yy); 

peak_time=(tp-1)*0.1;

overshoot=(yymax-1)*100




rr=1;


while y(rr)<1.0001
    
    
    rr=rr+1;
    
end 

rise_time=(rr-1)*0.1
    


s=(tf/2);

while y(s)>0.98 & y(s)<1.02;
    s=s-1;
end


settling_time=(s-1)*0.1



iae=0.1*sum(abs(e))


h=0.1;

for i= 1:length(t)
    
g(i)=0.01*i*e(i);

end

itae=sum(abs(g));

itae

  
