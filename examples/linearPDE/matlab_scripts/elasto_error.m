clear;
clc;

h = [1.0/3.0;
     0.2;
     0.125;
     0.1];
el2 = [5.431942e-07;
       1.157627e-07;
       2.747701e-08;
       1.393060e-08];
eh1 = [1.029172e-05;
       3.695816e-06;
       1.416088e-06;
       8.995649e-07];

h = log(h);
el2 = log(el2);
eh1 = log(eh1);

figure(1)
hold on
plot(h,el2,"-*",LineWidth=2);
plot(h,eh1,"-*","LineWidth",2);
legend("L2","H1","FontSize",14,"FontWeight","bold",Location="northwest");
title("Hex27","FontSize",14,"FontWeight","bold");
xlabel("ln(h)","FontSize",14,"FontWeight","bold");
ylabel("ln(e)","FontSize",14,"FontWeight","bold");
for i = 1:3
    text(0.5*(h(i)+h(i+1))-0.1,0.5*(el2(i)+el2(i+1))+0.1,num2str((el2(i)-el2(i+1))/(h(i)-h(i+1))),"FontSize",12,"FontWeight","bold");
    text(0.5*(h(i)+h(i+1))-0.1,0.5*(eh1(i)+eh1(i+1))+0.1,num2str((eh1(i)-eh1(i+1))/(h(i)-h(i+1))),"FontSize",12,"FontWeight","bold");
end
box on

el2 = [1.444874e-06;
       3.260426e-07;
       8.161793e-08;
       4.207758e-08];
eh1 = [3.664555e-05;
       1.444191e-05;
       5.860315e-06;
       3.786877e-06];

el2 = log(el2);
eh1 = log(eh1);

figure(2)
hold on
plot(h,el2,"-*",LineWidth=2);
plot(h,eh1,"-*","LineWidth",2);
legend("L2","H1","FontSize",14,"FontWeight","bold",Location="northwest");
title("Tet10","FontSize",14,"FontWeight","bold");
xlabel("ln(h)","FontSize",14,"FontWeight","bold");
ylabel("ln(e)","FontSize",14,"FontWeight","bold");
for i = 1:3
    text(0.5*(h(i)+h(i+1))-0.1,0.5*(el2(i)+el2(i+1))+0.1,num2str((el2(i)-el2(i+1))/(h(i)-h(i+1))),"FontSize",12,"FontWeight","bold");
    text(0.5*(h(i)+h(i+1))-0.1,0.5*(eh1(i)+eh1(i+1))+0.1,num2str((eh1(i)-eh1(i+1))/(h(i)-h(i+1))),"FontSize",12,"FontWeight","bold");
end
box on

h = [1.0/3.0;
     0.2;
     0.1;
     0.05];
el2 = [4.765113e-6;
       1.675566e-6;
       4.201252e-7;
       1.052733e-7];
eh1 = [7.195680e-5;
       4.242117e-5;
       2.109046e-5;
       1.053089e-5];

h = log(h);
el2 = log(el2);
eh1 = log(eh1);

figure(3)
hold on
plot(h,el2,"-*",LineWidth=2);
plot(h,eh1,"-*","LineWidth",2);
legend("L2","H1","FontSize",14,"FontWeight","bold",Location="northwest");
title("Hex8","FontSize",14,"FontWeight","bold");
xlabel("ln(h)","FontSize",14,"FontWeight","bold");
ylabel("ln(e)","FontSize",14,"FontWeight","bold");
for i = 1:3
    text(0.5*(h(i)+h(i+1))-0.1,0.5*(el2(i)+el2(i+1))+0.1,num2str((el2(i)-el2(i+1))/(h(i)-h(i+1))),"FontSize",12,"FontWeight","bold");
    text(0.5*(h(i)+h(i+1))-0.1,0.5*(eh1(i)+eh1(i+1))+0.1,num2str((eh1(i)-eh1(i+1))/(h(i)-h(i+1))),"FontSize",12,"FontWeight","bold");
end
box on

el2 = [1.182127e-5;
       4.814144e-6;
       1.301888e-6;
       3.370888e-7];
eh1 = [1.344155e-4;
       8.834327e-5;
       4.584205e-5;
       2.311859e-5];

el2 = log(el2);
eh1 = log(eh1);

figure(4)
hold on
plot(h,el2,"-*",LineWidth=2);
plot(h,eh1,"-*","LineWidth",2);
legend("L2","H1","FontSize",14,"FontWeight","bold",Location="northwest");
title("Tet4","FontSize",14,"FontWeight","bold");
xlabel("ln(h)","FontSize",14,"FontWeight","bold");
ylabel("ln(e)","FontSize",14,"FontWeight","bold");
for i = 1:3
    text(0.5*(h(i)+h(i+1))-0.1,0.5*(el2(i)+el2(i+1))+0.1,num2str((el2(i)-el2(i+1))/(h(i)-h(i+1))),"FontSize",12,"FontWeight","bold");
    text(0.5*(h(i)+h(i+1))-0.1,0.5*(eh1(i)+eh1(i+1))+0.1,num2str((eh1(i)-eh1(i+1))/(h(i)-h(i+1))),"FontSize",12,"FontWeight","bold");
end
box on