function [residual, g1, g2, g3] = CP_mon_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(40, 1);
T58 = exp(y(31))*(1-params(9))*exp(y(43))*exp(y(23))/exp(y(26));
T77 = params(3)*exp(y(59))*(1-params(8))*(exp(y(60))-exp(y(34)))+exp(y(59))*params(3)*params(8)*exp(y(61))*exp(y(64));
T117 = exp(y(41))*exp(y(40))*(1-y(57))/(exp(y(10))*(1-y(17)));
T153 = exp(y(53))*exp(y(26))^(1-params(9));
T162 = (exp(y(54))*exp(y(45))*exp(y(3)))^params(9);
T167 = exp(y(24))*params(9)*exp(y(43))/exp(y(45));
T174 = exp(y(3))*exp(y(54))*params(29)*(exp(y(45))-1)^params(7);
T187 = exp(y(24))*params(9)*exp(y(43))/exp(y(3))+exp(y(54))*(exp(y(30))-exp(y(51)));
T200 = (y(52)+params(31))/(params(31)+y(13))-1;
T206 = params(10)/2*T200^2;
T214 = exp(y(59))*params(3)*params(10)*((params(31)+y(67))/(y(52)+params(31))-1);
T215 = ((params(31)+y(67))/(y(52)+params(31)))^2;
T216 = T214*T215;
T244 = params(12)*exp(y(55))^params(11);
T250 = exp(y(16))^(params(11)*(-params(13)));
T254 = T244*T250*exp(y(11));
T258 = params(12)*exp(y(16))^(params(13)*(1-params(12)));
T260 = exp(y(55))^(params(12)-1);
T263 = (1-T258*T260)/(1-params(12));
T281 = exp(y(59))*params(3)*params(12)*exp(y(55))^(params(13)*(-params(11)));
T284 = exp(y(68))^params(11);
T288 = T281*T284*exp(y(65));
T296 = exp(y(59))*params(3)*params(12)*exp(y(55))^(params(13)*(1-params(11)));
T298 = exp(y(68))^(params(11)-1);
T302 = T296*T298*exp(y(66));
T310 = params(11)*exp(y(55))*exp(y(48))/exp(y(49))/(params(11)-1);
T326 = exp(y(12))^params(16);
T330 = 1/params(3)*exp(y(55))^params(14);
T332 = exp(y(46))/(params(11)/(params(11)-1));
T334 = T332^params(15);
T335 = T330*T334;
T337 = T335^(1-params(16));
T341 = T326*T337*exp(y(22));
lhs =exp(y(31));
rhs =(exp(y(28))-params(5)*exp(y(4)))^(-params(4))-params(5)*params(3)*(exp(y(58))-exp(y(28))*params(5))^(-params(4));
residual(1)= lhs-rhs;
lhs =params(3)*exp(y(59));
rhs =1/exp(y(34));
residual(2)= lhs-rhs;
lhs =exp(y(32));
rhs =exp(y(31))/exp(y(7));
residual(3)= lhs-rhs;
lhs =params(28)*exp(y(26))^params(6);
rhs =T58;
residual(4)= lhs-rhs;
lhs =exp(y(38));
rhs =T77;
residual(5)= lhs-rhs;
lhs =exp(y(39));
rhs =1-params(8)+exp(y(59))*params(3)*params(8)*exp(y(62))*exp(y(63));
residual(6)= lhs-rhs;
lhs =exp(y(40));
rhs =exp(y(39))/((params(27)-exp(y(38)))*(1-y(57)));
residual(7)= lhs-rhs;
lhs =exp(y(41));
rhs =exp(y(8))+exp(y(10))*(1-y(17))*(exp(y(33))-exp(y(8)));
residual(8)= lhs-rhs;
lhs =exp(y(42));
rhs =T117;
residual(9)= lhs-rhs;
lhs =exp(y(18));
rhs =exp(y(40))*exp(y(35));
residual(10)= lhs-rhs;
lhs =exp(y(35));
rhs =exp(y(37))+exp(y(36));
residual(11)= lhs-rhs;
lhs =exp(y(36));
rhs =exp(y(41))*params(8)*exp(y(9));
residual(12)= lhs-rhs;
lhs =exp(y(37));
rhs =(1-y(17))*params(26)*exp(y(1));
residual(13)= lhs-rhs;
lhs =exp(y(18));
rhs =exp(y(30))*exp(y(25));
residual(14)= lhs-rhs;
lhs =exp(y(24));
rhs =T153*T162;
residual(15)= lhs-rhs;
lhs =T167;
rhs =T174;
residual(16)= lhs-rhs;
lhs =exp(y(44));
rhs =exp(y(43))*(1-params(9))*exp(y(23))/exp(y(26));
residual(17)= lhs-rhs;
lhs =exp(y(33));
rhs =T187/exp(y(6));
residual(18)= lhs-rhs;
lhs =exp(y(30))-1;
rhs =(y(52)+params(31))*params(10)*T200/(params(31)+y(13))+T206-T216;
residual(19)= lhs-rhs;
lhs =exp(y(51));
rhs =params(30)+params(29)*(exp(y(45))-1)^(1+params(7))/(1+params(7));
residual(20)= lhs-rhs;
lhs =y(52);
rhs =exp(y(27))-exp(y(3))*exp(y(54))*exp(y(51));
residual(21)= lhs-rhs;
lhs =exp(y(25));
rhs =y(52)+exp(y(54))*exp(y(3));
residual(22)= lhs-rhs;
lhs =exp(y(24));
rhs =exp(y(23))*exp(y(47));
residual(23)= lhs-rhs;
lhs =exp(y(47));
rhs =T254+(1-params(12))*T263^((-params(11))/(1-params(12)));
residual(24)= lhs-rhs;
lhs =exp(y(46));
rhs =1/exp(y(43));
residual(25)= lhs-rhs;
lhs =exp(y(48));
rhs =exp(y(43))*exp(y(23))+T288;
residual(26)= lhs-rhs;
lhs =exp(y(49));
rhs =exp(y(23))+T302;
residual(27)= lhs-rhs;
lhs =exp(y(56));
rhs =T310;
residual(28)= lhs-rhs;
lhs =exp(y(55))^(1-params(11));
rhs =(1-params(12))*exp(y(56))^(1-params(11))+params(12)*exp(y(16))^(params(13)*(1-params(11)));
residual(29)= lhs-rhs;
lhs =exp(y(50));
rhs =exp(y(34))*exp(y(68));
residual(30)= lhs-rhs;
lhs =exp(y(50));
rhs =T341;
residual(31)= lhs-rhs;
lhs =exp(y(23));
rhs =exp(y(28))+exp(y(27))+exp(y(18))*y(57)*params(25)+exp(y(29))+(y(52)+params(31))*T206;
residual(32)= lhs-rhs;
lhs =y(29);
rhs =params(21)*y(5)+x(it_, 3);
residual(33)= lhs-rhs;
lhs =y(57);
rhs =y(17)*params(23)+x(it_, 4);
residual(34)= lhs-rhs;
lhs =y(53);
rhs =params(19)*y(14)+x(it_, 1);
residual(35)= lhs-rhs;
lhs =y(54);
rhs =params(17)*y(15)+x(it_, 2);
residual(36)= lhs-rhs;
lhs =y(22);
rhs =params(1)*y(2)+x(it_, 5);
residual(37)= lhs-rhs;
lhs =y(19);
rhs =y(23)-(steady_state(6));
residual(38)= lhs-rhs;
lhs =y(21);
rhs =y(28)-(steady_state(11));
residual(39)= lhs-rhs;
lhs =y(20);
rhs =y(27)-(steady_state(10));
residual(40)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(40, 73);

  %
  % Jacobian matrix
  %

T400 = (-(exp(y(43))*(1-params(9))*exp(y(23))/exp(y(26))));
T405 = (-(exp(y(24))*params(9)*exp(y(43))/exp(y(3))/exp(y(6))));
T409 = (-(T153*exp(y(54))*exp(y(45))*exp(y(3))*getPowerDeriv(exp(y(54))*exp(y(45))*exp(y(3)),params(9),1)));
T440 = getPowerDeriv(exp(y(28))-params(5)*exp(y(4)),(-params(4)),1);
T445 = getPowerDeriv(exp(y(58))-exp(y(28))*params(5),(-params(4)),1);
T531 = getPowerDeriv(T335,1-params(16),1);
T564 = params(10)/2*(-(y(52)+params(31)))/((params(31)+y(13))*(params(31)+y(13)))*2*T200;
T575 = params(10)/2*2*T200*1/(params(31)+y(13));
T613 = getPowerDeriv(T263,(-params(11))/(1-params(12)),1);
  g1(1,4)=(-((-(params(5)*exp(y(4))))*T440));
  g1(1,28)=(-(exp(y(28))*T440-params(5)*params(3)*(-(exp(y(28))*params(5)))*T445));
  g1(1,58)=params(5)*params(3)*exp(y(58))*T445;
  g1(1,31)=exp(y(31));
  g1(2,59)=params(3)*exp(y(59));
  g1(2,34)=(-((-exp(y(34)))/(exp(y(34))*exp(y(34)))));
  g1(3,7)=(-((-(exp(y(31))*exp(y(7))))/(exp(y(7))*exp(y(7)))));
  g1(3,31)=(-(exp(y(31))/exp(y(7))));
  g1(3,32)=exp(y(32));
  g1(4,23)=(-T58);
  g1(4,26)=params(28)*exp(y(26))*getPowerDeriv(exp(y(26)),params(6),1)-(-(exp(y(26))*exp(y(31))*(1-params(9))*exp(y(43))*exp(y(23))))/(exp(y(26))*exp(y(26)));
  g1(4,31)=(-T58);
  g1(4,43)=(-T58);
  g1(5,59)=(-T77);
  g1(5,60)=(-(params(3)*exp(y(59))*(1-params(8))*exp(y(60))));
  g1(5,34)=(-(params(3)*exp(y(59))*(1-params(8))*(-exp(y(34)))));
  g1(5,38)=exp(y(38));
  g1(5,61)=(-(exp(y(59))*params(3)*params(8)*exp(y(61))*exp(y(64))));
  g1(5,64)=(-(exp(y(59))*params(3)*params(8)*exp(y(61))*exp(y(64))));
  g1(6,59)=(-(exp(y(59))*params(3)*params(8)*exp(y(62))*exp(y(63))));
  g1(6,39)=exp(y(39));
  g1(6,62)=(-(exp(y(59))*params(3)*params(8)*exp(y(62))*exp(y(63))));
  g1(6,63)=(-(exp(y(59))*params(3)*params(8)*exp(y(62))*exp(y(63))));
  g1(7,38)=(-((-(exp(y(39))*(1-y(57))*(-exp(y(38)))))/((params(27)-exp(y(38)))*(1-y(57))*(params(27)-exp(y(38)))*(1-y(57)))));
  g1(7,39)=(-(exp(y(39))/((params(27)-exp(y(38)))*(1-y(57)))));
  g1(7,40)=exp(y(40));
  g1(7,57)=(-((-(exp(y(39))*(-(params(27)-exp(y(38))))))/((params(27)-exp(y(38)))*(1-y(57))*(params(27)-exp(y(38)))*(1-y(57)))));
  g1(8,33)=(-(exp(y(10))*(1-y(17))*exp(y(33))));
  g1(8,8)=(-(exp(y(8))+exp(y(10))*(1-y(17))*(-exp(y(8)))));
  g1(8,10)=(-(exp(y(10))*(1-y(17))*(exp(y(33))-exp(y(8)))));
  g1(8,41)=exp(y(41));
  g1(8,17)=(-((exp(y(33))-exp(y(8)))*(-exp(y(10)))));
  g1(9,10)=(-(exp(y(41))*(-(exp(y(10))*(1-y(17))*exp(y(40))*(1-y(57))))/(exp(y(10))*(1-y(17))*exp(y(10))*(1-y(17)))));
  g1(9,40)=(-T117);
  g1(9,41)=(-T117);
  g1(9,42)=exp(y(42));
  g1(9,17)=(-(exp(y(41))*(-(exp(y(40))*(1-y(57))*(-exp(y(10)))))/(exp(y(10))*(1-y(17))*exp(y(10))*(1-y(17)))));
  g1(9,57)=(-(exp(y(41))*(-exp(y(40)))/(exp(y(10))*(1-y(17)))));
  g1(10,18)=exp(y(18));
  g1(10,35)=(-(exp(y(40))*exp(y(35))));
  g1(10,40)=(-(exp(y(40))*exp(y(35))));
  g1(11,35)=exp(y(35));
  g1(11,36)=(-exp(y(36)));
  g1(11,37)=(-exp(y(37)));
  g1(12,9)=(-(exp(y(41))*params(8)*exp(y(9))));
  g1(12,36)=exp(y(36));
  g1(12,41)=(-(exp(y(41))*params(8)*exp(y(9))));
  g1(13,1)=(-((1-y(17))*params(26)*exp(y(1))));
  g1(13,37)=exp(y(37));
  g1(13,17)=params(26)*exp(y(1));
  g1(14,18)=exp(y(18));
  g1(14,25)=(-(exp(y(30))*exp(y(25))));
  g1(14,30)=(-(exp(y(30))*exp(y(25))));
  g1(15,24)=exp(y(24));
  g1(15,3)=T409;
  g1(15,26)=(-(T162*exp(y(53))*exp(y(26))*getPowerDeriv(exp(y(26)),1-params(9),1)));
  g1(15,45)=T409;
  g1(15,53)=(-(T153*T162));
  g1(15,54)=T409;
  g1(16,24)=T167;
  g1(16,3)=(-T174);
  g1(16,43)=T167;
  g1(16,45)=(-(exp(y(45))*exp(y(24))*params(9)*exp(y(43))))/(exp(y(45))*exp(y(45)))-exp(y(3))*exp(y(54))*params(29)*exp(y(45))*getPowerDeriv(exp(y(45))-1,params(7),1);
  g1(16,54)=(-T174);
  g1(17,23)=T400;
  g1(17,26)=(-(exp(y(43))*(-(exp(y(26))*(1-params(9))*exp(y(23))))/(exp(y(26))*exp(y(26)))));
  g1(17,43)=T400;
  g1(17,44)=exp(y(44));
  g1(18,24)=T405;
  g1(18,3)=(-((-(exp(y(3))*exp(y(24))*params(9)*exp(y(43))))/(exp(y(3))*exp(y(3)))/exp(y(6))));
  g1(18,6)=(-((-(T187*exp(y(6))))/(exp(y(6))*exp(y(6)))));
  g1(18,30)=(-(exp(y(30))*exp(y(54))/exp(y(6))));
  g1(18,33)=exp(y(33));
  g1(18,43)=T405;
  g1(18,51)=(-(exp(y(54))*(-exp(y(51)))/exp(y(6))));
  g1(18,54)=(-(exp(y(54))*(exp(y(30))-exp(y(51)))/exp(y(6))));
  g1(19,30)=exp(y(30));
  g1(19,59)=T216;
  g1(19,13)=(-(((params(31)+y(13))*(y(52)+params(31))*params(10)*(-(y(52)+params(31)))/((params(31)+y(13))*(params(31)+y(13)))-(y(52)+params(31))*params(10)*T200)/((params(31)+y(13))*(params(31)+y(13)))+T564));
  g1(19,52)=(-((params(10)*T200+(y(52)+params(31))*params(10)*1/(params(31)+y(13)))/(params(31)+y(13))+T575-(T215*exp(y(59))*params(3)*params(10)*(-(params(31)+y(67)))/((y(52)+params(31))*(y(52)+params(31)))+T214*(-(params(31)+y(67)))/((y(52)+params(31))*(y(52)+params(31)))*2*(params(31)+y(67))/(y(52)+params(31)))));
  g1(19,67)=T215*exp(y(59))*params(3)*params(10)*1/(y(52)+params(31))+T214*2*(params(31)+y(67))/(y(52)+params(31))*1/(y(52)+params(31));
  g1(20,45)=(-(params(29)*exp(y(45))*getPowerDeriv(exp(y(45))-1,1+params(7),1)/(1+params(7))));
  g1(20,51)=exp(y(51));
  g1(21,3)=exp(y(3))*exp(y(54))*exp(y(51));
  g1(21,27)=(-exp(y(27)));
  g1(21,51)=exp(y(3))*exp(y(54))*exp(y(51));
  g1(21,52)=1;
  g1(21,54)=exp(y(3))*exp(y(54))*exp(y(51));
  g1(22,3)=(-(exp(y(54))*exp(y(3))));
  g1(22,25)=exp(y(25));
  g1(22,52)=(-1);
  g1(22,54)=(-(exp(y(54))*exp(y(3))));
  g1(23,23)=(-(exp(y(23))*exp(y(47))));
  g1(23,24)=exp(y(24));
  g1(23,47)=(-(exp(y(23))*exp(y(47))));
  g1(24,11)=(-T254);
  g1(24,47)=exp(y(47));
  g1(24,16)=(-(exp(y(11))*T244*exp(y(16))*getPowerDeriv(exp(y(16)),params(11)*(-params(13)),1)+(1-params(12))*(-(T260*params(12)*exp(y(16))*getPowerDeriv(exp(y(16)),params(13)*(1-params(12)),1)))/(1-params(12))*T613));
  g1(24,55)=(-(exp(y(11))*T250*params(12)*exp(y(55))*getPowerDeriv(exp(y(55)),params(11),1)+(1-params(12))*T613*(-(T258*exp(y(55))*getPowerDeriv(exp(y(55)),params(12)-1,1)))/(1-params(12))));
  g1(25,43)=(-((-exp(y(43)))/(exp(y(43))*exp(y(43)))));
  g1(25,46)=exp(y(46));
  g1(26,23)=(-(exp(y(43))*exp(y(23))));
  g1(26,59)=(-T288);
  g1(26,43)=(-(exp(y(43))*exp(y(23))));
  g1(26,48)=exp(y(48));
  g1(26,65)=(-T288);
  g1(26,55)=(-(exp(y(65))*T284*exp(y(59))*params(3)*params(12)*exp(y(55))*getPowerDeriv(exp(y(55)),params(13)*(-params(11)),1)));
  g1(26,68)=(-(exp(y(65))*T281*exp(y(68))*getPowerDeriv(exp(y(68)),params(11),1)));
  g1(27,23)=(-exp(y(23)));
  g1(27,59)=(-T302);
  g1(27,49)=exp(y(49));
  g1(27,66)=(-T302);
  g1(27,55)=(-(exp(y(66))*T298*exp(y(59))*params(3)*params(12)*exp(y(55))*getPowerDeriv(exp(y(55)),params(13)*(1-params(11)),1)));
  g1(27,68)=(-(exp(y(66))*T296*exp(y(68))*getPowerDeriv(exp(y(68)),params(11)-1,1)));
  g1(28,48)=(-T310);
  g1(28,49)=(-(params(11)*(-(exp(y(49))*exp(y(55))*exp(y(48))))/(exp(y(49))*exp(y(49)))/(params(11)-1)));
  g1(28,55)=(-T310);
  g1(28,56)=exp(y(56));
  g1(29,16)=(-(params(12)*exp(y(16))*getPowerDeriv(exp(y(16)),params(13)*(1-params(11)),1)));
  g1(29,55)=exp(y(55))*getPowerDeriv(exp(y(55)),1-params(11),1);
  g1(29,56)=(-((1-params(12))*exp(y(56))*getPowerDeriv(exp(y(56)),1-params(11),1)));
  g1(30,34)=(-(exp(y(34))*exp(y(68))));
  g1(30,50)=exp(y(50));
  g1(30,68)=(-(exp(y(34))*exp(y(68))));
  g1(31,22)=(-T341);
  g1(31,46)=(-(exp(y(22))*T326*T330*T332*getPowerDeriv(T332,params(15),1)*T531));
  g1(31,12)=(-(exp(y(22))*T337*exp(y(12))*getPowerDeriv(exp(y(12)),params(16),1)));
  g1(31,50)=exp(y(50));
  g1(31,55)=(-(exp(y(22))*T326*T531*T334*1/params(3)*exp(y(55))*getPowerDeriv(exp(y(55)),params(14),1)));
  g1(32,18)=(-(exp(y(18))*y(57)*params(25)));
  g1(32,23)=exp(y(23));
  g1(32,27)=(-exp(y(27)));
  g1(32,28)=(-exp(y(28)));
  g1(32,29)=(-exp(y(29)));
  g1(32,13)=(-((y(52)+params(31))*T564));
  g1(32,52)=(-(T206+(y(52)+params(31))*T575));
  g1(32,57)=(-(exp(y(18))*params(25)));
  g1(33,5)=(-params(21));
  g1(33,29)=1;
  g1(33,71)=(-1);
  g1(34,17)=(-params(23));
  g1(34,57)=1;
  g1(34,72)=(-1);
  g1(35,14)=(-params(19));
  g1(35,53)=1;
  g1(35,69)=(-1);
  g1(36,15)=(-params(17));
  g1(36,54)=1;
  g1(36,70)=(-1);
  g1(37,2)=(-params(1));
  g1(37,22)=1;
  g1(37,73)=(-1);
  g1(38,19)=1;
  g1(38,23)=(-1);
  g1(39,21)=1;
  g1(39,28)=(-1);
  g1(40,20)=1;
  g1(40,27)=(-1);

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],40,5329);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],40,389017);
end
end
end
end
