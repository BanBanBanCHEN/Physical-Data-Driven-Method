function [y]=plot_error(path)
%%
    mkdir(['plot/',path]);
    Sim_U=struct();
    fid=fopen('plot/out.dat');
    while ~feof(fid)
        line=fgetl(fid);
        S=regexp(line,':','split');
        name=S{1};
        Ua=regexp(S{2},'\s+','split');
        u=[];
        for i=1:length(Ua)
            u=[u;str2double(Ua{i})];
        end
        Sim_U.(name)=u;
    end
    fclose(fid);
    %%
    fid=fopen('object.txt','r');
    Object=struct();
    while ~feof(fid)
        line=fgetl(fid);
        W=regexp(line,':','split');
        S=regexp(W{2},' ','split'); 
        u=[];
        for i =1:length(S)
            u=[u;str2double(S{i})];
        end
        Object.(W{1})=u;
    end
    fclose(fid);
    step=struct();
    step.M4008=[1,2,3,4,5,6,7,8,10,11]';
    step.M4067=[4,5,6,7,8,10,11]';
    step.M4068=[4,5,6,7,8,11]';
    A=fieldnames(Sim_U);
    for i=1:length(fieldnames(Sim_U))
        name=A{i};
        u=Sim_U.(name);
        obj=Object.(name);
        figure(i);
        plot(step.(name),u,'r','linewidth',1.5,'marker','o','markersize',10,'markerfacecolor','r',...
            'markeredgecolor','r','DisplayName','u');
        hold on
        plot(step.(name),obj,'k','linewidth',1.5,'marker','o','markersize',10,'markerfacecolor','k',...
            'markeredgecolor','k','DisplayName','Object');
        set(gca,'fontname','Times New Roman','fontsize',25,'fontweight','bold','linewidth',1.5);
        title(['\fontname{Times New Roman}',name,'\fontname{ËÎÌå}µÄÎ»ÒÆ'],'fontsize',30);
        xlabel('steps','fontsize',25);
        ylabel('u/mm','fontsize',25);
        a=max(max(u),max(obj));
        ylim([0,a]);
        set(gcf,'position',[0,0,1000,1000]);
        axis square;
        h=legend('show');
        h.Location='northeastoutside';
        print(i,'-dpng','-r300',['plot/',path,'/',name,'.png']);
    end
    close all;
    y=1;
end