x = 0:0.005:1;%Size
y1 = x.*(1-x);%logistics
y2 = x.^(2/3)-x;%von Ber...
y3 = -x.*log(x);%Gompertz
plot(x,y1, 'r:', x, y2, x, y3,'k-.')
%legend('logistic,von Bertalanffy,Gompertz',1)
ylim([0 0.4]);
xlabel('S')
ylabel("Growth functions");
title("Plot growth functions","FontSize",12);
