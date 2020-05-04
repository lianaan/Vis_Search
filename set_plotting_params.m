%function [guttera, marginsa, colors_sz, colors_sz_shade ] = set_plotting_params()

guttera=[0.14 0.06];
marginsa=[0.14 0.14 0.04 0.04]; %left right bottom top

%fontsz = 10;
fontsz = 13;
sz_dot = 1.2;

col_gray = [168 168 168]/255;
color_gray=[210 210 210]/255;
colors_sz = [217 30 0 ;  51 164 0;  242 135 65; 130 44 169]./255; % red, green, orange, purple
colors_sz_shade =(colors_sz+1.5*repmat(color_gray,4,1))/2;
colors_sz_shade(colors_sz_shade>0.99)= 0.99;

redd = [0.9047    0.1918    0.1988];
bluee = [0.2941    0.5447    0.7494];
redd_shade=(redd+2*[1 1 1])/3;
bluee_shade=(bluee+2*[1 1 1])/3;

text_labels = {'Perception', 'Memory'};

%end