



opt.objfunction = 'osy';
opt = basic_parameters(opt);
pop = lhsamp_model(100, opt);
[pop_obj, pop_cons] = evaluate_pop(opt, pop);
pop_cv = evaluateCV(opt, pop_cons);

dlmwrite(opt.varfilename, pop, 'delimiter',' ','precision','%.10f','-append');
dlmwrite(opt.objfilename, pop_obj, 'delimiter',' ','precision','%.10f','-append');
dlmwrite(opt.cvfilename, pop_cv, 'delimiter', ' ','precision','%.10f','-append');

