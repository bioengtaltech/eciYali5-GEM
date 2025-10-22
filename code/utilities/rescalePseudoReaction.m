function model = rescalePseudoReaction(model,metName,f)
  % rescalePseudoReaction
  %   Rescales a specific pseudoreaction by a given factor
  %
  %   model      (struct) the yeast GEM
  %   metName    (str) name of the component to rescale (e.g. "protein")
  %   f          (float) fraction to use for rescaling
  %
  %   model    (struct) the (rescaled) yeast GEM
  %
  %   Usage: model = rescalePseudoReaction(model,metName,f)
  %

rxnName = [metName ' pseudoreaction'];
rxnPos  = strcmp(model.rxnNames,rxnName);
for i = 1:length(model.mets)
    S_ir   = model.S(i,rxnPos);
    isProd = strcmp(model.metNames{i},[metName ' [cytoplasm]']);
    if S_ir ~= 0 && ~isProd
        S_ir = full(S_ir);
        model.S(i,rxnPos) = f*S_ir;
    end
end

end