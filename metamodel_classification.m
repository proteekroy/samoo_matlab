% MIT License
% 
% Copyright (c) 2018 Proteek Roy
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

% Contact:
% Proteek Chandan Roy
% Email: proteek_buet@yahoo.com


function metamodel_classification()
        
    clear all;
    rng('shuffle');
%     test_function = ['zdt1   ';'zdt2   ';'zdt3   ';'zdt6   ';...
%                      'osy    ';'tnk    ';'srn    ';'bnh    ';...
%                      'welded ';'dtlz2  ';'c2dtlz2'; 'carside';...
%                      'dtlz5  ';'dtlz4  ';'dtlz7  ';'zdt4   ';];
    test_function = ['zdt1   ';'zdt2   ';'zdt3   ';'zdt4   ';'zdt6   ';...
                     'dtlz1  ';'dtlz2  ';'dtlz3  ';'dtlz4  ';...
                     'dtlz5  ';'dtlz7  ';'c2dtlz2';'osy    ';...
                     'tnk    ';'welded ';'srn    ';'bnh    ';...
                     'dtlz2  '; 'carside';...
                     ];
    
%-------------------MAIN LOOP----------------------------------------------

    methods = [11,12,21,22,31,32,41,42,5,6];
    
    
    
    for methodology = 9:1:9 %for all methodologies
        
        for func_no = 12:1:12% for all test functions  
            
            for r = 1:1%for all runs
                
                opt.r = r;%run number
                opt.objfunction = strtrim(test_function(func_no,:));%remove whitespaces
                disp(['Problem Name: ' upper(opt.objfunction)]);
                opt.func_no = func_no;%function number
                opt.methodology = methods(methodology); %methodology
                
                %--------------------SET BASIC PARAMETERS------------------
                opt = basic_parameters(opt);
                
                %-----------------------OPTIMIZE---------------------------
                opt = parallel_metamodel_main_loop(opt);
                clear opt;
                
            end  
        end
    end
    delete(gcp('nocreate'))
end


%------------------------------END OF -FILE--------------------------------

