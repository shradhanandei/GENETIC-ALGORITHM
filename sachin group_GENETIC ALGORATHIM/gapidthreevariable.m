%%%%%exp(-0.2s) /(s+1)^2 

clc ;

clear all;


popsize=10;

Nt=12;


no_of_variable=3;


pop=round(rand(popsize,Nt));


aa=pop(:,1:4);

bb=pop(:,5:8); 

cc=pop(:,9:12);


d1=bi2de(aa,'left-msb');
d2=bi2de(bb,'left-msb');
d3=bi2de(cc,'left-msb');


%%%initialization %%%%%%%%%%%
  
    
ku=10.5;
tu=2.0333;
 
kp=0.6*ku;

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

w1=1;
w2=0;

nmut=2;


h=0.1;   

t = 0:h:16; 


tf=160

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:10
       

x1(i,1)=xl1+((xh1-xl1)/(2^4-1))*d1(i);

x2(i,1)=xl2+((xh2-xl2)/(2^4-1))*d2(i);

x3(i,1)=xl3+((xh3-xl3)/(2^4-1))*d3(i);


end


d=[x1 x2 x3 ] ;

for j=1:10
    
    
    kp=d(j,1);
    ki=d(j,2);
    kd=d(j,3);
    
    

% % % 
% % % h=0.1;   
% % % 
% % % % step size
% % % 
% % % t = 0:h:16; 
% % % tf=160


y = zeros(1,length(t)); 
u = zeros(1,length(t));
e = zeros(1,length(t)); 
x = zeros(1,length(t));
input = 1;
y(1) = 0;
x(1)=0 ;
e(1)=input - y(1);

u(1)=kp*e(1)+ki*sum(e);

F_xy = @(x) -2*x; 

for i = 1:2
  %y and u taken as time input and x as output
    
  
   k_1 = F_xy(x(i));
   k_2 = F_xy(x(i)+0.5*h*k_1);
   k_3 = F_xy((x(i)+0.5*h*k_2));
   k_4 = F_xy((x(i)+k_3*h));

   
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
    y(i+1)=y(i)+h*x(i+1);
    
   e(i+1)=1-y(i+1);
   
 
   ed=e(i+1)-e(i);
   
   
 
   u(i+1)=kp*e(i+1)+ki*sum(e)+kd*ed;
   

end


F_xy = @(u,y,x) u-y-2*x;


for i = 3:80
    
    k_1 = F_xy(u(i-2),y(i),x(i));
    
    k_2 = F_xy(u(i-2)+0.5*h,y(i)+0.5*h,x(i)+0.5*h*k_1);
    
    k_3 = F_xy((u(i-2)+0.5*h),(y(i)+0.5*h),(x(i)+0.5*h*k_2));
    
    k_4 = F_xy((u(i-2)+h ) ,(y(i)+h) ,(x(i)+k_3*h));
    
    
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
   y(i+1)=y(i)+h*x(i+1);
   
   e(i+1)=1-y(i+1);
  
   
   ed=e(i+1)-e(i);
   
  u(i+1)=kp*e(i+1)+ki*sum(e)+kd*ed;
   
   
end


u(81)=-20;

F_xy = @(u,y,x) u-y-2*x;


for i = 81:length(t)
    k_1 = F_xy(u(i-2),y(i),x(i));
    
    k_2 = F_xy(u(i-2)+0.5*h,y(i)+0.5*h,x(i)+0.5*h*k_1);
    
    k_3 = F_xy((u(i-2)+0.5*h),(y(i)+0.5*h),(x(i)+0.5*h*k_2));
    
    k_4 = F_xy((u(i-2)+h ) ,(y(i)+h) ,(x(i)+k_3*h));
    
    
   
   x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
   y(i+1)=y(i)+h*x(i+1);
   
   e(i+1)=1-y(i+1);
  
   
   ed=e(i+1)-e(i);
   
  
   u(i+1)=kp*e(i+1)+ki*sum(e)+kd*ed;
    
end

z=y(1:161);

figure(j+1)

plot(t,z) 


xlabel('Time t(s) ')

ylabel('Response y')

title(' exp(-0.2s) /(s+1)^2  ')

grid on 

%%%%iaeitae

iaeiae(j,1)=0.1*sum(abs(e));


h=0.1;

for i= 1:length(t)
    
g(i)=0.01*i*e(i);

end

itaeitae(j,1)=sum(abs(g));

 end

iae=w1*iaeiae+w2*itaeitae;



[iae ind]=sort(iae);

pop=pop(ind,:);


d=d(ind,:);


iga=0;
 
 
 while iga<50
     
    
    

iga=iga+1;

%%%%%%%%%%%%%%%%selection of mimimum five iae%%%%%%%%%%%%%%%%%%%%%%

  
  pop;

  p1=[pop(1,1:2) pop(2,3:4) pop(1,5:6) pop(2,7:8) pop(1,9:10) pop(2,11:12) ];        
 
  p2=[pop(2,1:2) pop(1,3:4) pop(2,5:6) pop(1,7:8) pop(2,9:10) pop(1,11:12) ];   
 
  p3=[pop(3,1:2) pop(4,3:4) pop(3,5:6) pop(4,7:8) pop(3,9:10) pop(4,11:12)  ];        
  
  p4=[pop(4,1:2) pop(3,3:4) pop(4,5:6) pop(3,7:8) pop(4,9:10) pop(3,11:12)  ];        
  
  p5=[pop(5,1:2) pop(1,3:4) pop(5,5:6) pop(1,7:8) pop(5,9:10) pop(1,11:12)  ];        
  
  
  
  
  
  m =[p1;p2;p3;p4;p5];
  
  pop(6:10,:)=m ;
  
  
aa=pop(:,1:4);

bb=pop(:,5:8); 

cc=pop(:,9:12);

d1=bi2de(aa,'left-msb');
d2=bi2de(bb,'left-msb');
d3=bi2de(cc,'left-msb');
  


for i=1:10
       

x1(i,1)=xl1+((xh1-xl1)/(2^4-1))*d1(i);

x2(i,1)=xl2+((xh2-xl2)/(2^4-1))*d2(i);

x3(i,1)=xl3+((xh3-xl3)/(2^4-1))*d3(i);




end

d=[x1 x2 x3 ] ;

  
 %%%%objective function 
 

for j=1:10
    
    
    kp=d(j,1);
    ki=d(j,2);
    kd=d(j,3);
    
   
y = zeros(1,length(t)); 
u = zeros(1,length(t));
e = zeros(1,length(t)); 
x = zeros(1,length(t));
input = 1;
y(1) = 0;
x(1)=0 ;
e(1)=input - y(1);
  
  
 u(1)=kp*e(1)+ki*sum(e);

F_xy = @(x) -2*x; 

for i = 1:2
    
    
  %y and u taken as time input and x as output
   
   k_1 = F_xy(x(i));
   k_2 = F_xy(x(i)+0.5*h*k_1);
   k_3 = F_xy((x(i)+0.5*h*k_2));
   k_4 = F_xy((x(i)+k_3*h));

    

   
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
    y(i+1)=y(i)+h*x(i+1);
    
   e(i+1)=1-y(i+1);
   
   ed=e(i+1)-e(i);
   
  
   u(i+1)=kp*e(i+1)+ki*sum(e)+kd*ed;
   

end


F_xy = @(u,y,x) u-y-2*x;


for i = 3:80
    
    k_1 = F_xy(u(i-2),y(i),x(i));
    
    k_2 = F_xy(u(i-2)+0.5*h,y(i)+0.5*h,x(i)+0.5*h*k_1);
    
    k_3 = F_xy((u(i-2)+0.5*h),(y(i)+0.5*h),(x(i)+0.5*h*k_2));
    
    k_4 = F_xy((u(i-2)+h ) ,(y(i)+h) ,(x(i)+k_3*h));
    
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
   y(i+1)=y(i)+h*x(i+1);
   
   e(i+1)=1-y(i+1);
  
   
   ed=e(i+1)-e(i);
   
  u(i+1)=kp*e(i+1)+ki*sum(e)+kd*ed;
   
   
end


u(81)=-20



F_xy = @(u,y,x) u-y-2*x;


for i = 81:length(t)
    
    k_1 = F_xy(u(i-2),y(i),x(i));
    
    k_2 = F_xy(u(i-2)+0.5*h,y(i)+0.5*h,x(i)+0.5*h*k_1);
    
    k_3 = F_xy((u(i-2)+0.5*h),(y(i)+0.5*h),(x(i)+0.5*h*k_2));
    
    k_4 = F_xy((u(i-2)+h ) ,(y(i)+h) ,(x(i)+k_3*h));
    
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
   y(i+1)=y(i)+h*x(i+1);
   
   e(i+1)=1-y(i+1);
  
   
   ed=e(i+1)-e(i);
   
  
   
   u(i+1)=kp*e(i+1)+ki*sum(e)+kd*ed;
    
end



z=y(1:161);

figure(j+1)
plot(t,z) 


xlabel('Time t(s) ')

ylabel('Response y')

title(' exp(-0.2s) /(s+1)^2  ')

grid on 



%%%%iaeitae

iaeiae(j,1)=0.1*sum(abs(e))


h=0.1;

for i= 1:length(t)
    
g(i)=0.01*i*e(i);

end

itaeitae(j,1)=sum(abs(g));



 end



iae=w1*iaeiae+w2*itaeitae

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




pop ;

aa=pop(:,1:4);

bb=pop(:,5:8); 

cc=pop(:,9:12);

d1=bi2de(aa,'left-msb');
d2=bi2de(bb,'left-msb');
d3=bi2de(cc,'left-msb');
  


for i=1:10
       

x1(i,1)=xl1+((xh1-xl1)/(2^4-1))*d1(i);

x2(i,1)=xl2+((xh2-xl2)/(2^4-1))*d2(i);

x3(i,1)=xl3+((xh3-xl3)/(2^4-1))*d3(i);




end




d=[x1 x2 x3 ] ;

  

 
for j=1:10
    
    
    kp=d(j,1);
    ki=d(j,2);
    kd=d(j,3);
    
    


h=0.1;   

% step size

t = 0:h:16; 
tf=160


y = zeros(1,length(t)); 
u = zeros(1,length(t));
e = zeros(1,length(t)); 
x = zeros(1,length(t));
input = 1;
y(1) = 0;
x(1)=0 ;
e(1)=input - y(1);

  
 u(1)=kp*e(1)+ki*sum(e);


F_xy = @(x) -2*x; 

for i = 1:2
  %y and u taken as time input and x as output
    
  
  
   k_1 = F_xy(x(i));
   k_2 = F_xy(x(i)+0.5*h*k_1);
   k_3 = F_xy((x(i)+0.5*h*k_2));
   k_4 = F_xy((x(i)+k_3*h));

    

   
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
    y(i+1)=y(i)+h*x(i+1);
    
   e(i+1)=1-y(i+1);
   
   
   
   ed=e(i+1)-e(i);
   
   
 
   u(i+1)=kp*e(i+1)+ki*sum(e)+kd*ed;
   


end


F_xy = @(u,y,x) u-y-2*x;


for i = 3:80
    
    k_1 = F_xy(u(i-2),y(i),x(i));
    
    k_2 = F_xy(u(i-2)+0.5*h,y(i)+0.5*h,x(i)+0.5*h*k_1);
    
    k_3 = F_xy((u(i-2)+0.5*h),(y(i)+0.5*h),(x(i)+0.5*h*k_2));
    
    k_4 = F_xy((u(i-2)+h ) ,(y(i)+h) ,(x(i)+k_3*h));
    
    
    
    
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
   y(i+1)=y(i)+h*x(i+1);
   
   e(i+1)=1-y(i+1);
  
   
   ed=e(i+1)-e(i);
   
  u(i+1)=kp*e(i+1)+ki*sum(e)+kd*ed;
   
   
end


u(81)=-20



F_xy = @(u,y,x) u-y-2*x;


for i = 81:length(t)
    
    k_1 = F_xy(u(i-2),y(i),x(i));
    
    k_2 = F_xy(u(i-2)+0.5*h,y(i)+0.5*h,x(i)+0.5*h*k_1);
    
    k_3 = F_xy((u(i-2)+0.5*h),(y(i)+0.5*h),(x(i)+0.5*h*k_2));
    
    k_4 = F_xy((u(i-2)+h ) ,(y(i)+h) ,(x(i)+k_3*h));
    
    
    
    
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
   y(i+1)=y(i)+h*x(i+1);
   
   e(i+1)=1-y(i+1);
  
   
   ed=e(i+1)-e(i);
   
  
   
   u(i+1)=kp*e(i+1)+ki*sum(e)+kd*ed;
    
end



z=y(1:161);

figure(j+1)
plot(t,z) 


xlabel('Time t(s) ')

ylabel('Response y')

title(' exp(-0.2s) /(s+1)^2  ')

grid on 



%%%%iaeitae

iaeiae(j,1)=0.1*sum(abs(e))


h=0.1;

for i= 1:length(t)
    
g(i)=0.01*i*e(i);

end

itaeitae(j,1)=sum(abs(g));



 end



iae=w1*iaeiae+w2*itaeitae

[iae ind]=sort(iae);

pop=pop(ind,:);

d=d(ind,:)


 end 



ind;

iga

d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
 for j=1
    kp=d(j,1);
    ki=d(j,2);
    kd=d(j,3);
    
    
% % % 
% % % 
% % % h=0.1;   
% % % 
% % % % step size
% % % 
% % % t = 0:h:16; 



y = zeros(1,length(t)); 
u = zeros(1,length(t));
e = zeros(1,length(t)); 
x = zeros(1,length(t));
input = 1;
y(1) = 0;
x(1)=0 ;
e(1)=input - y(1);

  
 u(1)=kp*e(1)+ki*sum(e);

F_xy = @(x) -2*x; 

for i = 1:2
    
  %y and u taken as time input and x as output
    
  
  
   k_1 = F_xy(x(i));
   k_2 = F_xy(x(i)+0.5*h*k_1);
   k_3 = F_xy((x(i)+0.5*h*k_2));
   k_4 = F_xy((x(i)+k_3*h));

    

   
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
    y(i+1)=y(i)+h*x(i+1);
    
   e(i+1)=1-y(i+1);
   
   
   
   ed=e(i+1)-e(i);
   
   
 
   u(i+1)=kp*e(i+1)+ki*sum(e)+kd*ed;
   


end


F_xy = @(u,y,x) u-y-2*x;


for i = 3:80
    
    k_1 = F_xy(u(i-2),y(i),x(i));
    
    k_2 = F_xy(u(i-2)+0.5*h,y(i)+0.5*h,x(i)+0.5*h*k_1);
    
    k_3 = F_xy((u(i-2)+0.5*h),(y(i)+0.5*h),(x(i)+0.5*h*k_2));
    
    k_4 = F_xy((u(i-2)+h ) ,(y(i)+h) ,(x(i)+k_3*h));
    
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
   y(i+1)=y(i)+h*x(i+1);
   
   e(i+1)=1-y(i+1);
  
   
   ed=e(i+1)-e(i);
   
  u(i+1)=kp*e(i+1)+ki*sum(e)+kd*ed;
   
   
end


u(81)=-20;



F_xy = @(u,y,x) u-y-2*x;


for i = 81:length(t)
    
    k_1 = F_xy(u(i-2),y(i),x(i));
    
    k_2 = F_xy(u(i-2)+0.5*h,y(i)+0.5*h,x(i)+0.5*h*k_1);
    
    k_3 = F_xy((u(i-2)+0.5*h),(y(i)+0.5*h),(x(i)+0.5*h*k_2));
    
    k_4 = F_xy((u(i-2)+h ) ,(y(i)+h) ,(x(i)+k_3*h));
     
    x(i+1) = x(i) + (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
    
   y(i+1)=y(i)+h*x(i+1);
   
   e(i+1)=1-y(i+1);
  
   
   ed=e(i+1)-e(i);
   
  
   u(i+1)=kp*e(i+1)+ki*sum(e)+kd*ed;
    
end

z=y(1:161);

figure(1)

plot(t,z) 

xlabel('Time t(s) ')

ylabel('Response y')

title(' exp(-0.2s) /(s+1)^2  ')

grid on 

 end


 
%%%%% time domain specification -over shoot,risetime ,settling time
%%%%% ,iae,itae


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

 