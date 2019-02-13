% Copyright [2016] [Proteek Chandan Roy]
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.

%
% Proteek Chandan Roy, 2016
% Contact: royprote@msu.edu, proteek_buet@yahoo.com


function [mutated_pop, nrealmut] = pol_mut(pop_crossover, pmut, nrealmut, eta_m, lowerBound, upperBound)



[N, nreal] = size( pop_crossover); % Population size & Number of variables
mutated_pop =  pop_crossover;% Child before mutation


for ind = 1:N
    for i = 1:nreal
        if rand <= pmut
            y = mutated_pop(ind,i);
            yl = lowerBound(i);
            yu = upperBound(i);
            delta1 = (y-yl) / (yu-yl);
            delta2 = (yu-y) / (yu-yl);
            random = rand;
            mut_pow = 1.0/(eta_m+1.0);
            if random <= 0.5
                xy = 1.0 - delta1;
                val = 2.0*random + (1.0 - 2.0*random) * xy^(eta_m+1.0);
                deltaq =  val^mut_pow - 1.0;
            else
                xy = 1.0 - delta2;
                val = 2.0*(1.0 - random) + 2.0*(random-0.5) * xy^(eta_m+1.0);
                deltaq = 1.0 - val^mut_pow;
            end
            y = y + deltaq*(yu - yl);
            if (y<yl), y = yl; end
            if (y>yu), y = yu; end
            mutated_pop(ind,i) = y;
            nrealmut = nrealmut+1;
        end
    end
end

%------------------------------END OF -FILE--------------------------------
