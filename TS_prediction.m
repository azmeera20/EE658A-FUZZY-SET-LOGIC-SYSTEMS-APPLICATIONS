clear;close all;clc;
% loading Mackey-Glass Chaotic Time Series
load mgdata.dat; % Inbuilt data
m = 9; % points taken for predicting one value
figure(1);
plot(mgdata(1:1000,1),mgdata(1:1000,2),'LineWidth',1.5);
xlabel('t');ylabel('x(t)');
title('Mackey-Glass Chaotic Time Series for {\tau} = 30');

input_mgdata = 0;
output_mgdata = 0;
for i = 1:(size(mgdata,1)-m) % creating (M-m) pairs of input and output
    input_mgdata(i,1:m) = mgdata(i:i+m-1,2)'; % taking only second column because first column has only indices to represent on x-axis.
    output_mgdata(i,1) = mgdata(i+m,2)'; % taking (m+1)th column value as output which we will get in next state.
end
%% Train 700 data used
train_data = [input_mgdata(1:700,:) output_mgdata(1:700,:)]; % according to paper training data has 700 samples.
test_data1 = input_mgdata(701:1000,:); % according to paper test data has remaining 300 samples.
N=[4 14];
for i=1:2
    test_output(i,:)=train_test(train_data,test_data1,N(i));
    disp('mean square error=');
    norm(test_output(i,:)-output_mgdata(701:1000))
end 

figure(3);
plot(1:300,output_mgdata(701:1000),'LineWidth',2.1);hold on;
plot(1:300,test_output(1,:),'LineWidth',1.7);hold on
plot(1:300,test_output(2,:),'LineWidth',1.5)
title('Chaotic Time Series Prediction when 700 train data used');
xlabel('x(t)','FontSize',12);
ylabel('\mu(x)','FontSize',12);
legend('Actual Value','Predicted Value N=4','Predicted Value N=14');
%% Train 200 data used
train_data = [input_mgdata(501:700,:) output_mgdata(501:700,:)];% according to paper training data has 200 samples.
test_data2 = [input_mgdata(1:150,:); input_mgdata(851:1000,:)];% according to paper test data has 300 samples.
for i=1:2
    test_output(i,:)=train_test(train_data,test_data2,N(i));
    disp('mean square error=')
    norm(test_output(i,:)-output_mgdata(701:1000))
end

figure(5);
actual_output=[output_mgdata(1:150);output_mgdata(851:1000)];
plot(1:300,actual_output,'LineWidth',2.1);hold on;
plot(1:300,test_output(1,:),'LineWidth',1.7);hold on;
plot(1:300,test_output(2,:),'LineWidth',1.5)
title('Chaotic Time Series Prediction when 200 train data used');
xlabel('x(t)','FontSize',12);
ylabel('\mu(x)','FontSize',12);
legend('Actual Value','Predicted Value at N=4','Predicted Value at N=14');
%% Sequence of steps followed here
function test_output=train_test(train_data,test_data,N)
    m=9;
    %% Step 1: Divide the input and output spaces into fuzzy regions.
    % Divide each domain interval into 2N+1 regions
   % N = 4;
    X = 0.1:0.01:2.3;
    [X_FuzzyReg, R] = FuzzyRegions(N, X);
    y_t = 0; % y_t Finds center value with membership value = 1
        for j = 1:size(X_FuzzyReg,2)
            [x_loc,y_loc] = find(X_FuzzyReg{:,j}==1);
            y_t(1,j) = X(1,y_loc);
            figure(N);
            plot(X,X_FuzzyReg{1,j},'Linewidth',1.5);
            hold on;
            title(sprintf('Membership Function for Chaotic Time Series Prediction N=%d,R=%d',N,R));
            xlabel('x(t)');ylabel('\mu(x)');ylim([0 1]);
        end
    %% Step 2: Generate Fuzzy Rules from Given Data Pairs.
    Degree_Value = 0;
    Rule_Value = 0;
    for k = 1:size(train_data,1)
        for l = 1:size(train_data,2)
            degree_data = zeros(1,R);
            train_data(k,l) = round(train_data(k,l)*100)/100;
            for J = 1:R
                temp1 = repmat(train_data(k,l),size(X));
                temp2 = temp1 - X;
                temp3 = abs(temp2);
                indx = find(temp3 < 0.001);
                if (indx)
                    degree_data(1,J) = X_FuzzyReg{1,J}(indx);
                end
            end
                [Degree_Value(k,l), Rule_Value(k,l)] = max(degree_data);
         end
    end
    %% Step 3: Assign a Degree to Each Rule.
    % For conflicting rule in rule base, assign the maximum degree for the
    Degree_Rule1 = prod(Degree_Value,2);
    [tmp, index] = unique(Rule_Value,'rows','stable'); % unique rules. (max among the similar will be taken)
    NewRule_Degree = Degree_Rule1(index); % Degree of unique rules obtained.
    [a,b,c] = unique(tmp(:,1:m),'rows','stable'); % value of m is 9 as given in paper (see before step 1).
    n = 1;
    % final_matrix = zeros(1,1);
    for i = 1:size(b,1)
        dup_rows_index = find(c==i); % Identify rows having same index
        if length(dup_rows_index) > 1 % Check no. of matching
            [u,v] = max(NewRule_Degree(dup_rows_index)); % Find row having maximum rule degree
            final_matrix(n,:) = tmp(dup_rows_index(v),:); % Keep that rule and put in final matrix
        else
            final_matrix(n,:) = tmp(dup_rows_index,:) ;
        end
        n = n + 1;
    end
    %% Step 4: Create a Combined Fuzzy Rule Base.
    fuzzy_rule_base = final_matrix;
    % Step 5: Determine a Mapping Based on the Combined Fuzzy Rule Base.
    Degree_test = 0;
    in_mf_prod = 0;
    y_bar = zeros(1,1);
    test_output = 0;
        for p = 1:size(test_data,1)
            sample = test_data(p,:);
            for q = 1:size(fuzzy_rule_base,1)
                for r = 1:size(test_data,2)
                    val = zeros(1,R);
                    train_data = round(sample(r)*100)/100;
                    for I = 1:R
                        temp1 = repmat(train_data,size(X));
                        temp2 = temp1 - X;
                        temp3 = abs(temp2);
                        indx = find(temp3 < 0.001);
                        if (indx)
                            val(1,I) = X_FuzzyReg{1,I}(indx);
                        end
                    end
                    Degree_test(1,r) = val(fuzzy_rule_base(q,r));
                end
            in_mf_prod(q,:) = prod(Degree_test);
            y_bar(q,:) = y_t(fuzzy_rule_base(q,10)) ;
            end
        test_output(p,1) = sum(in_mf_prod.*y_bar)/sum(in_mf_prod) ;
        end
end