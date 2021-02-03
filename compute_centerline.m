function compute_centerline(dir, M, structure, orden, decimate, optim, subj_name, hemi)
%compute_centerline uses function from hippovol
%(https://github.com/garikoitz/hippovol) to compute centerline of a region
% INPUTS: 
%      - dir: directory where output images will be saved
%      - M: volume obtained with MRIread
%      - structure: name of structure being analysed
%      - orden: order for Bezier function
%      - decimate: decimation in Bezier function
%      - optim: 1: use matlab's internal fmninunc, 0: use lbfsg in cluster 
%      - subj_name: name of subject being processed
%      - hemi: lh: left hemisphere, rh: right hemisphere
%

    [ShortBezier,ctrl_pts] = hip_fitBezierToCloud_vGari(M, orden, decimate, optim, subj_name);
    save([dir,'/',subj_name,'/',subj_name,'_',hemi,'_',structure,'_bezier_pts.txt'],'ShortBezier','-ascii','-double','-tabs')
    save([dir,'/',subj_name,'/',subj_name,'_',hemi,'_',structure,'_ctrl_pts.txt'],'ctrl_pts','-ascii','-double','-tabs')

    %p=patch(isosurface(M.vol),'visible','off'); isonormals(M.vol, p)
    
    stlwrite([dir,'/',subj_name,'/',subj_name,'_',hemi,'_',structure,'.stl'],isosurface(M.vol,64.64))
%     value = isovalue(M.vol);
%     sprintf('Isovalue: %.2f', value)
end
