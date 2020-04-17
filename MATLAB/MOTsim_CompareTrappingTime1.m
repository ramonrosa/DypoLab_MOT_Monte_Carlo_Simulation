function MOTsim_CompareTrappingTime1 ()

    [DT1,TTm1] = MOTsim_ViewResults1(0);
    [DT2,TTm2] = MOTsim_ViewResults1(0);


    %%

    plot(DT1,TTm1,'Color',[0 0 0],'Marker','o','MarkerFaceColor',[0 0 0.5]);
    hold on;
    plot(DT2,TTm2,'Color',[0 0 0],'Marker','s','MarkerFaceColor',[1 1 1]);
    hold off;
    set(gca,'XDir','Reverse');
    xlabel('Detuning (units of gamma)');
    ylabel('Trapping time (s)');
    legend('Configuration A','Configuration B');
    set(gca,'XTick',-150:25:0);
    grid on;