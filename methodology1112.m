function opt = methodology1112(opt)


    opt.switching_option = 1;
    opt = methodology11(opt);
    opt.switching_option = 2;
    opt = nsga2_main(opt);

end