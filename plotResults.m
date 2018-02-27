function []=plotResults(results)
	figure(3)
	time = results.time;
	symbols = {'r-o','b-*'};
    for wellId = 1:length(results.wellQs)
        subplot(3,1,1);
        qo = abs(results.wellQs{wellId}(:,1));
        plot(time,qo,symbols{wellId});
        title('oil rate');
        hold on;
        subplot(3,1,2);
        qw = abs(results.wellQs{wellId}(:,2));
        plot(time,qw,symbols{wellId});
        title('water rate');
        hold on;
        subplot(3,1,3);
        pbh = results.wellPbh{wellId}(:);
        plot(time,pbh,symbols{wellId});
        title('BHP');
        hold on;
    end
    for wellId = 1:length(results.wellQs)
        subplot(3,1,1);
        hold off;
        subplot(3,1,2);
        hold off;
        subplot(3,1,3);
        hold off;
    end
end
