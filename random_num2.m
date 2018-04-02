clear all;
close all;

%% CMRG和combRecursive源随机数

%参数
para = [32,1403580,-810728,-209,527612,-1370589,-22853];
 
%10000个随机数
generate_length = 10000;
 
%设置初始值
a1 = zeros(1,generate_length + 3);
a2 = zeros(1,generate_length + 3);
y = 0;
c = zeros(1,generate_length + 3);
d1 = zeros(1,generate_length);
 
%给定递推式开始值
a1(1) = 1; a1(2) = 2; a1(3) = 3;
a2(1) = 4; a2(2) = 5; a2(3) = 6;
 
for i = 4 : generate_length + 3
    a1(i) = mod(para(2) .* a1(i-2) + para(3) .* a1(i-3), 2 .^ para(1) + para(4));
    a2(i) = mod(para(5) .* a2(i-2) + para(6) .* a2(i-3), 2 .^ para(1) + para(7));
    y = mod(a1(i) - a2(i), 2 .^ para(1) + para(4));
    c(i) = y ./ (2 .^ para(1) + para(4));
end
 
%去掉生成数列中前三个值
for i = 1 : generate_length
    d1(i) = c(i + 3);%生成CMRG
end
 
rng('shuffle', 'combRecursive');% 修改generator源为'combRecursive'
d2 = rand(1,generate_length); %生成10000个随机数

figure(1)
hist(d1)
figure(2)
hist(d2)

%% 1、均匀性检验

k = 100;%k = 100个子区间
N1 = hist(d1,k);
N2 = hist(d2,k);
xd1 = zeros(1,k);
xd2 = zeros(1,k);

%等分成K个子区间
for i = 1 : k
    xd1(i) = (N1(i) - generate_length ./ k) .^ 2 ./ (generate_length ./ k);
    xd2(i) = (N2(i) - generate_length ./ k) .^ 2 ./ (generate_length ./ k);
end
 
kf_1_1 = chi2cdf(sum(xd1),k-1)
kf_1_2 = chi2cdf(sum(xd2),k-1)
 
%在原假设条件下，统计量接近K-1自由度的卡方分布
if kf_1_1 < 0.95
    kf_1_1_result = 'pass'
else
    kf_1_1_result = 'fail'
end
 
if kf_1_2 < 0.95
    kf_1_2_result = 'pass'
else
    kf_1_2_result = 'fail'
end

%% 2、序列检验

k = 64;
Hist1 = zeros(k);
Hist2 = zeros(k);

%构建2维卡方统计量
for i = 1 : 1 : length(d1) / 2
   x1_1(i) = d1(i);
   x1_2(i) = d1(i+1);
   x2_1(i) = d2(i);
   x2_2(i) = d2(i+1);
end
 
%将单位正方形分成K*K个相等的小正方形
for i = 1 : length(x1_1)
    j1 = ceil(x1_1(i)*k);
    j2 = ceil(x1_2(i)*k);
    j3 = ceil(x2_1(i)*k);
    j4 = ceil(x2_2(i)*k);
    Hist1(j1,j2) = Hist1(j1,j2) + 1;
    Hist2(j3,j4) = Hist2(j3,j2) + 1;
end

figure(3)
bar3(Hist1)
figure(4)
bar3(Hist2)
 
 
for i = 1:k
    for j = 1:k
        h1(i,j) = ((Hist1(i,j) - (generate_length / 2) / (k ^ 2)) .^ 2);
        h2(i,j) = ((Hist2(i,j) - (generate_length / 2) / (k ^ 2)) .^ 2);
    end
end
 
kf_2_1 = sum(sum(h1)) .* (k ^ 2) / (generate_length / 2)
kf_2_2 = sum(sum(h2)) .* (k ^ 2) / (generate_length / 2)
vpa(kf_2_1,5)%精度小数点后三位
vpa(kf_2_2,5)
 
if kf_2_1 < 4016.5
    kf_2_1_result = 'pass'
else
    kf_2_1_result = 'fail'
end
 
if kf_2_2 < 4016.5
    kf_2_2_result = 'pass'
else
    kf_2_2_result = 'fail'
end

%% 3、游程检验
 
a_ij = [4529.4	9044.9	13568	18091	22615	27892
        9044.9	18097	27139	36187	45234	55789
        13568   27139	40721	54281	67852	83685
        18091   36187	54281	72414	90470	111580
        22615   45234	67852	90470	113262	139476
        27892   55789	83685	111580	139476	172860];
 
b_i = [1/6,5/24,11/120,19/720,29/5040,1/840];
 
w = [1,2,3,4,5,6];
 
r_ij_1 = zeros(1,6);
r_ij_2 = zeros(1,6);
temp1 = 1;
temp2 = 1;
 
for i = 1 : length(d1)-1
    if d1(i+1)>d1(i)
        temp1 = temp1+1;
    else
        if temp1<6
            r_ij_1(temp1) = r_ij_1(temp1) + 1;
        else
            r_ij_1(6) = r_ij_1(6) + 1;
        end
        temp1 = 1;
    end
    
    if d2(i+1)>d2(i)
        temp2=temp2+1;
    else
        if temp2<6
            r_ij_2(temp2) = r_ij_2(temp2) + 1;
        else
            r_ij_2(6) = r_ij_2(6) + 1;
        end
        temp2=1;
    end
end
 
sum(r_ij_1.*w)
sum(r_ij_2.*w)
 
for i=1:6
    for j=1:6
        s1(i,j)=a_ij(i,j)*(r_ij_1(i)-generate_length*b_i(i))*(r_ij_1(j)-generate_length*b_i(j));
        s2(i,j)=a_ij(i,j)*(r_ij_2(i)-generate_length*b_i(i))*(r_ij_2(j)-generate_length*b_i(j));
    end
end
 
kf_3_1 = chi2cdf(sum(sum(s1)) / generate_length, 6)
kf_3_2 = chi2cdf(sum(sum(s2)) / generate_length, 6)
 
if kf_3_1 < 0.95
    kf_3_1_result = 'pass'
else
    kf_3_1_result = 'fail'
end
 
if kf_3_2 < 0.95
    kf_3_2_result = 'pass'
else
    kf_3_2_result = 'fail'
end

%% 4、相关性检验 

j=3;
h=(generate_length-1)/j-1;
for i=0:h
    r1(i+1)=d1(1+i*j)*d1(1+(i+1)*j);
    r2(i+1)=d2(1+i*j)*d2(1+(i+1)*j);
end
 
R1=sum(r1) * 12 / (h + 1) - 3;
v=(13 * h + 7) / ((h + 1) .^ 2);
A1=abs(R1 / sqrt(v))
 
normp1=normcdf(0,1,A1)
 
R2=sum(r2) * 12 / (h + 1) - 3;
v=(13 * h + 7) / ((h + 1) .^ 2);
A2=abs(R2 / sqrt(v))
 
normp2=normcdf(0,1,A2)
 
if normp1 < 0.95
    kf_4_1_result = 'pass'
else
    kf_4_1_result = 'fail'
end
 
if normp2 < 0.95
    kf_4_2_result = 'pass'
else
    kf_4_2_result = 'fail'
end

