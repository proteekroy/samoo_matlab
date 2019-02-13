function [net] = train_neural_network(x, t)


net = network;
net.numInputs = 1;
net.numLayers = 2;%3;
net.biasConnect = [1; 1];%[1;1;1];
net.inputConnect = [1; 0];%[1;0;0];
net.layerConnect = [0 0;1 0];%[0 0 0;1 0 0;0 1 0];%column denotes the input
net.outputConnect = [0 1];%[0 0 1];


net.adaptFcn = 'adaptwb';
net.divideFcn = 'dividerand'; %Set the divide function to dividerand (divide training data randomly).

net.performFcn = 'mse';
net.trainFcn = 'trainlm'; % set training function to trainlm (Levenberg-Marquardt backpropagation) 

net.plotFcns = {'plotperform', 'plottrainstate', 'ploterrhist', 'plotconfusion', 'plotroc'};

hidden_layer_size1 = 10;
hidden_layer_size2 = 10;

%set Layer1
net.layers{1}.name = 'Layer 1';
net.layers{1}.dimensions = hidden_layer_size1;
net.layers{1}.initFcn = 'initnw';
net.layers{1}.transferFcn = 'poslin'; %'tansig';

%set Layer2
net.layers{1}.name = 'Layer 2';
net.layers{1}.dimensions = hidden_layer_size2;%hidden_layer_size;
net.layers{1}.initFcn = 'initnw';
net.layers{1}.transferFcn = 'poslin'; %'tansig';

%set Layer3
% net.layers{2}.name = 'Layer 3';
% net.layers{2}.dimensions = 1; %hidden_layer_size;
% net.layers{2}.initFcn = 'initnw';
% net.layers{2}.transferFcn = 'purelin';%'poslin'; %'tansig';

net.trainParam.epochs = 20;
net.trainParam.showWindow = 0; %default is 1)
net = train(net,x', t'); %training

%y = net(x); %prediction

end