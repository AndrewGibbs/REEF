% Create figure window and components
function FFslider(theta_test,Emb_obj)
fig = uifigure('Position',[100 100 350 275]);
title('inc field slider');

sld = uislider(fig,...
    'Position',[100 75 120 3],...
    'ValueChangedFcn',@(sld,event) updateGauge(event));
sld.MajorTickLabels = {'0','pi/2','pi','3pi/2','2pi'};
sld.MajorTicks = [0 pi/2 pi 3*pi/2 2*pi];
sld.Limits = [0 2*pi];

    function updateGauge(event)
        Evals = Emb_obj.getFarField(theta_test,event.Value);
        plot(theta_test,real(Evals),theta_test,imag(Evals));
        ylim([-10*Emb_obj.kwave 10*Emb_obj.kwave]);
        legend('real','imaginary');
        title('Far-field pattern');
        xlim([0 2*pi]);
        xticks([0 pi/2 pi 3*pi/2 2*pi]);
        xticklabels({'0','pi/2','pi','3pi/2','2pi'});
        shg;
    end

end