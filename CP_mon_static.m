function [residual, g1, g2, g3] = CP_mon_static(y, x, params)
%
% Status : Computes static model for Dynare
%
% Inputs : 
%   y         [M_.endo_nbr by 1] double    vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1] double     vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1] double   vector of parameter values in declaration order
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the static model equations 
%                                          in order of declaration of the equations.
%                                          Dynare may prepend or append auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by M_.endo_nbr] double    Jacobian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g2        [M_.endo_nbr by (M_.endo_nbr)^2] double   Hessian matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%   g3        [M_.endo_nbr by (M_.endo_nbr)^3] double   Third derivatives matrix of the static model equations;
%                                                       columns: variables in declaration order
%                                                       rows: equations in order of declaration
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

residual = zeros( 40, 1);

%
% Model equations
%

T17 = (exp(y(11))-exp(y(11))*params(5))^(-params(4));
T46 = exp(y(14))*(1-params(9))*exp(y(26))*exp(y(6))/exp(y(9));
T63 = params(3)*exp(y(15))*(1-params(8))*(exp(y(16))-exp(y(17)))+exp(y(21))*exp(y(15))*params(3)*params(8)*exp(y(25));
T117 = exp(y(36))*exp(y(9))^(1-params(9));
T124 = (exp(y(8))*exp(y(37))*exp(y(28)))^params(9);
T129 = exp(y(7))*params(9)*exp(y(26))/exp(y(28));
T136 = exp(y(8))*exp(y(37))*params(29)*(exp(y(28))-1)^params(7);
T149 = exp(y(7))*params(9)*exp(y(26))/exp(y(8))+exp(y(37))*(exp(y(13))-exp(y(34)));
T183 = exp(y(38));
T185 = T183^params(11);
T190 = T183^(params(11)*(-params(13)));
T192 = exp(y(30))*params(12)*T185*T190;
T196 = params(12)*T183^(params(13)*(1-params(12)));
T198 = T183^(params(12)-1);
T219 = exp(y(15))*params(3)*params(12)*T183^(params(13)*(-params(11)));
T221 = exp(y(31))*T185*T219;
T228 = T183^(params(13)*(1-params(11)));
T231 = T183^(params(11)-1);
T233 = exp(y(32))*exp(y(15))*params(3)*params(12)*T228*T231;
T241 = params(11)*T183*exp(y(31))/exp(y(32))/(params(11)-1);
T254 = exp(y(33))^params(16);
T258 = 1/params(3)*T183^params(14);
T260 = exp(y(29))/(params(11)/(params(11)-1));
T262 = T260^params(15);
T263 = T258*T262;
T265 = T263^(1-params(16));
T269 = T254*T265*exp(y(5));
lhs =exp(y(14));
rhs =T17-T17*params(5)*params(3);
residual(1)= lhs-rhs;
lhs =params(3)*exp(y(15));
rhs =1/exp(y(17));
residual(2)= lhs-rhs;
lhs =exp(y(15));
rhs =1;
residual(3)= lhs-rhs;
lhs =params(28)*exp(y(9))^params(6);
rhs =T46;
residual(4)= lhs-rhs;
lhs =exp(y(21));
rhs =T63;
residual(5)= lhs-rhs;
lhs =exp(y(22));
rhs =1-params(8)+exp(y(15))*params(3)*params(8)*exp(y(22))*exp(y(24));
residual(6)= lhs-rhs;
lhs =exp(y(23));
rhs =exp(y(22))/((params(27)-exp(y(21)))*(1-y(40)));
residual(7)= lhs-rhs;
lhs =exp(y(24));
rhs =exp(y(17))+(exp(y(16))-exp(y(17)))*exp(y(23))*(1-y(40));
residual(8)= lhs-rhs;
lhs =exp(y(25));
rhs =exp(y(24));
residual(9)= lhs-rhs;
lhs =exp(y(1));
rhs =exp(y(23))*exp(y(18));
residual(10)= lhs-rhs;
lhs =exp(y(18));
rhs =exp(y(20))+exp(y(19));
residual(11)= lhs-rhs;
lhs =exp(y(19));
rhs =exp(y(24))*params(8)*exp(y(18));
residual(12)= lhs-rhs;
lhs =exp(y(20));
rhs =(1-y(40))*exp(y(1))*params(26);
residual(13)= lhs-rhs;
lhs =exp(y(1));
rhs =exp(y(13))*exp(y(8));
residual(14)= lhs-rhs;
lhs =exp(y(7));
rhs =T117*T124;
residual(15)= lhs-rhs;
lhs =T129;
rhs =T136;
residual(16)= lhs-rhs;
lhs =exp(y(27));
rhs =exp(y(26))*(1-params(9))*exp(y(6))/exp(y(9));
residual(17)= lhs-rhs;
lhs =exp(y(16));
rhs =T149/exp(y(13));
residual(18)= lhs-rhs;
residual(19) = exp(y(13))-1;
lhs =exp(y(34));
rhs =params(30)+params(29)*(exp(y(28))-1)^(1+params(7))/(1+params(7));
residual(20)= lhs-rhs;
lhs =y(35);
rhs =exp(y(10))-exp(y(8))*exp(y(37))*exp(y(34));
residual(21)= lhs-rhs;
lhs =exp(y(8));
rhs =y(35)+exp(y(8))*exp(y(37));
residual(22)= lhs-rhs;
lhs =exp(y(7));
rhs =exp(y(6))*exp(y(30));
residual(23)= lhs-rhs;
lhs =exp(y(30));
rhs =T192+(1-params(12))*((1-T196*T198)/(1-params(12)))^((-params(11))/(1-params(12)));
residual(24)= lhs-rhs;
lhs =exp(y(29));
rhs =1/exp(y(26));
residual(25)= lhs-rhs;
lhs =exp(y(31));
rhs =exp(y(26))*exp(y(6))+T221;
residual(26)= lhs-rhs;
lhs =exp(y(32));
rhs =exp(y(6))+T233;
residual(27)= lhs-rhs;
lhs =exp(y(39));
rhs =T241;
residual(28)= lhs-rhs;
lhs =T183^(1-params(11));
rhs =(1-params(12))*exp(y(39))^(1-params(11))+params(12)*T228;
residual(29)= lhs-rhs;
lhs =exp(y(33));
rhs =exp(y(17))*T183;
residual(30)= lhs-rhs;
lhs =exp(y(33));
rhs =T269;
residual(31)= lhs-rhs;
lhs =exp(y(6));
rhs =exp(y(11))+exp(y(10))+exp(y(1))*y(40)*params(25)+exp(y(12));
residual(32)= lhs-rhs;
lhs =y(12);
rhs =y(12)*params(21)+x(3);
residual(33)= lhs-rhs;
lhs =y(40);
rhs =y(40)*params(23)+x(4);
residual(34)= lhs-rhs;
lhs =y(36);
rhs =y(36)*params(19)+x(1);
residual(35)= lhs-rhs;
lhs =y(37);
rhs =y(37)*params(17)+x(2);
residual(36)= lhs-rhs;
lhs =y(5);
rhs =y(5)*params(1)+x(5);
residual(37)= lhs-rhs;
lhs =y(2);
rhs =y(6)-(y(6));
residual(38)= lhs-rhs;
lhs =y(4);
rhs =y(11)-(y(11));
residual(39)= lhs-rhs;
lhs =y(3);
rhs =y(10)-(y(10));
residual(40)= lhs-rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
if nargout >= 2,
  g1 = zeros(40, 40);

  %
  % Jacobian matrix
  %

T322 = (-(exp(y(26))*(1-params(9))*exp(y(6))/exp(y(9))));
T327 = (-(exp(y(7))*params(9)*exp(y(26))/exp(y(8))/exp(y(13))));
T332 = (-(T117*exp(y(8))*exp(y(37))*exp(y(28))*getPowerDeriv(exp(y(8))*exp(y(37))*exp(y(28)),params(9),1)));
T362 = (exp(y(11))-exp(y(11))*params(5))*getPowerDeriv(exp(y(11))-exp(y(11))*params(5),(-params(4)),1);
T431 = getPowerDeriv(T263,1-params(16),1);
T463 = T183*getPowerDeriv(T183,params(11),1);
T495 = T183*getPowerDeriv(T183,params(13)*(1-params(11)),1);
  g1(1,11)=(-(T362-params(5)*params(3)*T362));
  g1(1,14)=exp(y(14));
  g1(2,15)=params(3)*exp(y(15));
  g1(2,17)=(-((-exp(y(17)))/(exp(y(17))*exp(y(17)))));
  g1(3,15)=exp(y(15));
  g1(4,6)=(-T46);
  g1(4,9)=params(28)*exp(y(9))*getPowerDeriv(exp(y(9)),params(6),1)-(-(exp(y(9))*exp(y(14))*(1-params(9))*exp(y(26))*exp(y(6))))/(exp(y(9))*exp(y(9)));
  g1(4,14)=(-T46);
  g1(4,26)=(-T46);
  g1(5,15)=(-T63);
  g1(5,16)=(-(params(3)*exp(y(15))*(1-params(8))*exp(y(16))));
  g1(5,17)=(-(params(3)*exp(y(15))*(1-params(8))*(-exp(y(17)))));
  g1(5,21)=exp(y(21))-exp(y(21))*exp(y(15))*params(3)*params(8)*exp(y(25));
  g1(5,25)=(-(exp(y(21))*exp(y(15))*params(3)*params(8)*exp(y(25))));
  g1(6,15)=(-(exp(y(15))*params(3)*params(8)*exp(y(22))*exp(y(24))));
  g1(6,22)=exp(y(22))-exp(y(15))*params(3)*params(8)*exp(y(22))*exp(y(24));
  g1(6,24)=(-(exp(y(15))*params(3)*params(8)*exp(y(22))*exp(y(24))));
  g1(7,21)=(-((-(exp(y(22))*(1-y(40))*(-exp(y(21)))))/((params(27)-exp(y(21)))*(1-y(40))*(params(27)-exp(y(21)))*(1-y(40)))));
  g1(7,22)=(-(exp(y(22))/((params(27)-exp(y(21)))*(1-y(40)))));
  g1(7,23)=exp(y(23));
  g1(7,40)=(-((-(exp(y(22))*(-(params(27)-exp(y(21))))))/((params(27)-exp(y(21)))*(1-y(40))*(params(27)-exp(y(21)))*(1-y(40)))));
  g1(8,16)=(-(exp(y(16))*exp(y(23))*(1-y(40))));
  g1(8,17)=(-(exp(y(17))+exp(y(23))*(1-y(40))*(-exp(y(17)))));
  g1(8,23)=(-((exp(y(16))-exp(y(17)))*exp(y(23))*(1-y(40))));
  g1(8,24)=exp(y(24));
  g1(8,40)=(-((exp(y(16))-exp(y(17)))*(-exp(y(23)))));
  g1(9,24)=(-exp(y(24)));
  g1(9,25)=exp(y(25));
  g1(10,1)=exp(y(1));
  g1(10,18)=(-(exp(y(23))*exp(y(18))));
  g1(10,23)=(-(exp(y(23))*exp(y(18))));
  g1(11,18)=exp(y(18));
  g1(11,19)=(-exp(y(19)));
  g1(11,20)=(-exp(y(20)));
  g1(12,18)=(-(exp(y(24))*params(8)*exp(y(18))));
  g1(12,19)=exp(y(19));
  g1(12,24)=(-(exp(y(24))*params(8)*exp(y(18))));
  g1(13,1)=(-((1-y(40))*exp(y(1))*params(26)));
  g1(13,20)=exp(y(20));
  g1(13,40)=exp(y(1))*params(26);
  g1(14,1)=exp(y(1));
  g1(14,8)=(-(exp(y(13))*exp(y(8))));
  g1(14,13)=(-(exp(y(13))*exp(y(8))));
  g1(15,7)=exp(y(7));
  g1(15,8)=T332;
  g1(15,9)=(-(T124*exp(y(36))*exp(y(9))*getPowerDeriv(exp(y(9)),1-params(9),1)));
  g1(15,28)=T332;
  g1(15,36)=(-(T117*T124));
  g1(15,37)=T332;
  g1(16,7)=T129;
  g1(16,8)=(-T136);
  g1(16,26)=T129;
  g1(16,28)=(-(exp(y(28))*exp(y(7))*params(9)*exp(y(26))))/(exp(y(28))*exp(y(28)))-exp(y(8))*exp(y(37))*params(29)*exp(y(28))*getPowerDeriv(exp(y(28))-1,params(7),1);
  g1(16,37)=(-T136);
  g1(17,6)=T322;
  g1(17,9)=(-(exp(y(26))*(-(exp(y(9))*(1-params(9))*exp(y(6))))/(exp(y(9))*exp(y(9)))));
  g1(17,26)=T322;
  g1(17,27)=exp(y(27));
  g1(18,7)=T327;
  g1(18,8)=(-((-(exp(y(8))*exp(y(7))*params(9)*exp(y(26))))/(exp(y(8))*exp(y(8)))/exp(y(13))));
  g1(18,13)=(-((exp(y(13))*exp(y(13))*exp(y(37))-exp(y(13))*T149)/(exp(y(13))*exp(y(13)))));
  g1(18,16)=exp(y(16));
  g1(18,26)=T327;
  g1(18,34)=(-(exp(y(37))*(-exp(y(34)))/exp(y(13))));
  g1(18,37)=(-(exp(y(37))*(exp(y(13))-exp(y(34)))/exp(y(13))));
  g1(19,13)=exp(y(13));
  g1(20,28)=(-(params(29)*exp(y(28))*getPowerDeriv(exp(y(28))-1,1+params(7),1)/(1+params(7))));
  g1(20,34)=exp(y(34));
  g1(21,8)=exp(y(8))*exp(y(37))*exp(y(34));
  g1(21,10)=(-exp(y(10)));
  g1(21,34)=exp(y(8))*exp(y(37))*exp(y(34));
  g1(21,35)=1;
  g1(21,37)=exp(y(8))*exp(y(37))*exp(y(34));
  g1(22,8)=exp(y(8))-exp(y(8))*exp(y(37));
  g1(22,35)=(-1);
  g1(22,37)=(-(exp(y(8))*exp(y(37))));
  g1(23,6)=(-(exp(y(6))*exp(y(30))));
  g1(23,7)=exp(y(7));
  g1(23,30)=(-(exp(y(6))*exp(y(30))));
  g1(24,30)=exp(y(30))-T192;
  g1(24,38)=(-(exp(y(30))*(T190*params(12)*T463+params(12)*T185*T183*getPowerDeriv(T183,params(11)*(-params(13)),1))+(1-params(12))*(-(T198*params(12)*T183*getPowerDeriv(T183,params(13)*(1-params(12)),1)+T196*T183*getPowerDeriv(T183,params(12)-1,1)))/(1-params(12))*getPowerDeriv((1-T196*T198)/(1-params(12)),(-params(11))/(1-params(12)),1)));
  g1(25,26)=(-((-exp(y(26)))/(exp(y(26))*exp(y(26)))));
  g1(25,29)=exp(y(29));
  g1(26,6)=(-(exp(y(26))*exp(y(6))));
  g1(26,15)=(-T221);
  g1(26,26)=(-(exp(y(26))*exp(y(6))));
  g1(26,31)=exp(y(31))-T221;
  g1(26,38)=(-(exp(y(31))*(T219*T463+T185*exp(y(15))*params(3)*params(12)*T183*getPowerDeriv(T183,params(13)*(-params(11)),1))));
  g1(27,6)=(-exp(y(6)));
  g1(27,15)=(-T233);
  g1(27,32)=exp(y(32))-T233;
  g1(27,38)=(-(exp(y(32))*(T231*exp(y(15))*params(3)*params(12)*T495+exp(y(15))*params(3)*params(12)*T228*T183*getPowerDeriv(T183,params(11)-1,1))));
  g1(28,31)=(-T241);
  g1(28,32)=(-(params(11)*(-(exp(y(32))*T183*exp(y(31))))/(exp(y(32))*exp(y(32)))/(params(11)-1)));
  g1(28,38)=(-T241);
  g1(28,39)=exp(y(39));
  g1(29,38)=T183*getPowerDeriv(T183,1-params(11),1)-params(12)*T495;
  g1(29,39)=(-((1-params(12))*exp(y(39))*getPowerDeriv(exp(y(39)),1-params(11),1)));
  g1(30,17)=(-(exp(y(17))*T183));
  g1(30,33)=exp(y(33));
  g1(30,38)=(-(exp(y(17))*T183));
  g1(31,5)=(-T269);
  g1(31,29)=(-(exp(y(5))*T254*T258*T260*getPowerDeriv(T260,params(15),1)*T431));
  g1(31,33)=exp(y(33))-exp(y(5))*T265*exp(y(33))*getPowerDeriv(exp(y(33)),params(16),1);
  g1(31,38)=(-(exp(y(5))*T254*T431*T262*1/params(3)*T183*getPowerDeriv(T183,params(14),1)));
  g1(32,1)=(-(exp(y(1))*y(40)*params(25)));
  g1(32,6)=exp(y(6));
  g1(32,10)=(-exp(y(10)));
  g1(32,11)=(-exp(y(11)));
  g1(32,12)=(-exp(y(12)));
  g1(32,40)=(-(exp(y(1))*params(25)));
  g1(33,12)=1-params(21);
  g1(34,40)=1-params(23);
  g1(35,36)=1-params(19);
  g1(36,37)=1-params(17);
  g1(37,5)=1-params(1);
  g1(38,2)=1;
  g1(39,4)=1;
  g1(40,3)=1;
  if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
  end
if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],40,1600);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],40,64000);
end
end
end
end
