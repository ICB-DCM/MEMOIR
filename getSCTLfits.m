function [ parameters_SCTL ] = getSCTLfits( Model,Data,n_starts)

iD = 1;
for iData = 1:length(Data)
    
    parameters.name = arrayfun(@(x) ['phi_{' '}'],Model.exp{1}.ind_phi,'UniformOutput',false);
    parameters.min = -5*ones(size(parameters.name));
    parameters.max = -parameters.min;
    parameters.number = length(parameters.max);
    

    
    optionsMultistart = PestoOptions();
    optionsMultistart.n_starts = n_starts;
    optionsMultistart.comp_type = 'sequential';
    optionsMultistart.mode = 'visual';
    optionsMultistart.obj_type = 'log-posterior';
    
    if(isfield(Data{iData},'SCTL'))
        if(~isfield(Data{iData}.SCTL,'T'))
            Data{iData}.SCTL.T = zeros(0,1,size(Data{iData}.SCTL.Y,3));
        end
        tmp = arrayfun(@(x) any(~isnan(Data{iData}.SCTL.Y(:,:,x)),2),1:size(Data{iData}.SCTL.Y,3),'UniformOutput',false);
        Data{iData}.SCTL.ind_y = [tmp{:}];
        tmp = arrayfun(@(x) any(~isnan(Data{iData}.SCTL.T(:,:,x)),2),1:size(Data{iData}.SCTL.T,3),'UniformOutput',false);
        Data{iData}.SCTL.ind_t = [tmp{:}];
        for icl = 1:size(Data{iData}.SCTL.Y,3)
            objective_function = @(phi) objective_phi_wrapper(phi,Model,Data,iData,icl);
            
            parameters_SCTL{iD,icl} = getMultiStarts(parameters, objective_function, optionsMultistart);
        end
        
       iD = iD+1; 
    end
end


end