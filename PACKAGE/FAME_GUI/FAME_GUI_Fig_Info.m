h.fig_info  = figure(2); clf

set(h.fig_info, 'position',    [0 0 500 800],...
           'color',       bcolor,...
           'menubar',     'none',...
           'toolbar',     'none',...
           'name',        'FAME Information',...
           'numbertitle', 'off');
% Lattice information panel    
%% Set panels
h.info_pnl_lattice = uipanel('parent',         h.fig_info,...
                             'title',           'Lattice',...
                             'backgroundcolor', bcolor,...
                             'ForegroundColor', [0 0 0],...
                             'fontsize',        8,...
                             'position',        [0.02 0.45 0.96 0.5]);
h.info_pnl_mesh = uipanel('parent',         h.fig_info,...
                             'title',           'Mesh',...
                             'backgroundcolor', bcolor,...
                             'ForegroundColor', [0 0 0],...
                             'fontsize',        8,...
                             'position',        [0.02 0.23 0.96 0.2]);
h.info_pnl_material = uipanel('parent',         h.fig_info,...
                             'title',           'Material',...
                             'backgroundcolor', bcolor,...
                             'ForegroundColor', [0 0 0],...
                             'fontsize',        8,...
                             'position',        [0.02 0.01 0.96 0.2]);                         
%% Lattice panel                         
h.info_text_lattice_type_name = uicontrol('parent',  h.info_pnl_lattice,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'Lattice Type:',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.01 0.85 0.31 0.14],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on'); 
h.info_text_lattice_type_value = uicontrol('parent',  h.info_pnl_lattice,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'simple_cubic',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.32 0.85 0.65 0.14],...
                                     'HorizontalAlignment', 'center',...
                                     'visible',             'on');                                  
h.info_text_lattice_constant = uicontrol('parent',  h.info_pnl_lattice,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'Lattice Constant:',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.01 0.71 0.31 0.14],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on'); 
h.info_text_lattice_constant_name(7) = uicontrol('parent',  h.info_pnl_lattice,...
                                 'style',               'text',...
                                 'units',               'normalized',...
                                 'string',              'Permutation',...
                                 'backgroundcolor',     bcolor,...
                                 'ForegroundColor',     [0 0 0],...
                                 'fontsize',            8,...
                                 'position',            [0.32 0.71 0.31 0.14],...
                                 'HorizontalAlignment', 'right',...
                                 'visible',             'on'); 
h.info_text_lattice_constant_value(7) = uicontrol('parent',  h.info_pnl_lattice,...
                                 'style',               'text',...
                                 'units',               'normalized',...
                                 'string',              '[1 2 3]',...
                                 'backgroundcolor',     bcolor,...
                                 'ForegroundColor',     [0 0 0],...
                                 'fontsize',            8,...
                                 'position',            [0.64 0.71 0.31 0.14],...
                                 'HorizontalAlignment', 'center',...
                                 'visible',             'on');                                        
cc      = {'a:' 'b:' 'c:' '�\:' '�]:' '�^:'};
val     = [ 1 1 1 90 90 90];
for i = 1:6
    h.info_text_lattice_constant_name(i) = uicontrol('parent',  h.info_pnl_lattice,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              cc{i},...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.05+0.15*(i-1) 0.57 0.05 0.14],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on'); 
    h.info_text_lattice_constant_value(i) = uicontrol('parent',  h.info_pnl_lattice,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              num2str(val(i)),...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.1+0.15*(i-1) 0.57 0.1 0.14],...
                                     'HorizontalAlignment', 'center',...
                                     'visible',             'on');                                      
end     
                                 
h.info_text_lattice_vector = uicontrol('parent',  h.info_pnl_lattice,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'Lattice Vectors:',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.01 0.43 0.31 0.14],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on'); 
cc      = {'a1:' 'a2:' 'a3:'};
val     = [ 1 0 0;
            0 1 0;
            0 0 1];
for i = 1:3
    h.info_text_lattice_vector_name(i) = uicontrol('parent',  h.info_pnl_lattice,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              cc{i},...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.05+0.315*(i-1) 0.29 0.07 0.14],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on'); 
    h.info_text_lattice_vector_value(i) = uicontrol('parent',  h.info_pnl_lattice,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              ['[',num2str(val(1,i)),',',num2str(val(2,i)),',',num2str(val(3,i)),']'],...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.12+0.315*(i-1) 0.29 0.2 0.14],...
                                     'HorizontalAlignment', 'center',...
                                     'visible',             'on');                                      
end                                     
h.info_text_BZ_path = uicontrol('parent',  h.info_pnl_lattice,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'Bri. Zone Path:',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.01 0.15 0.31 0.14],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on'); 
h.info_text_BZ_path_value = uicontrol('parent',  h.info_pnl_lattice,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'GXMRG',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.32 0.15 0.67 0.14],...
                                     'HorizontalAlignment', 'center',...
                                     'visible',             'on');
h.info_text_partition_num = uicontrol('parent',  h.info_pnl_lattice,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'Partition Number:',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0 0.01 0.32 0.14],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on');
h.info_text_partition_num_val = uicontrol('parent',  h.info_pnl_lattice,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              '10',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.34 0.01 0.63 0.14],...
                                     'HorizontalAlignment', 'center',...
                                         'visible',             'on');
%% Mesh Panel
h.info_text_edge_len = uicontrol('parent',  h.info_pnl_mesh,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'Edge Length:',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0 0.835 0.32  0.165],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on');
cc = {'Lx:' 'Ly:' 'Lz:'};
for i = 1:3
    h.info_text_edge_len_name(i) = uicontrol('parent',  h.info_pnl_mesh,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              cc{i},...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.05+0.315*(i-1) 0.685 0.07 0.165],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on'); 
    h.info_text_edge_len_value(i) = uicontrol('parent',  h.info_pnl_mesh,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              '1',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.12+0.315*(i-1) 0.685 0.2 0.165],...
                                     'HorizontalAlignment', 'center',...
                                     'visible',             'on');    
end
h.info_text_grid_num = uicontrol('parent',  h.info_pnl_mesh,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'Grid Numbers:',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0 0.505 0.32  0.165],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on');
cc = {'Nx:' 'Ny:' 'Nz:'};
for i = 1:3
    h.info_text_grid_num_name(i) = uicontrol('parent',  h.info_pnl_mesh,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              cc{i},...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.05+0.315*(i-1) 0.34 0.07 0.165],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on'); 
    h.info_text_grid_num_value(i) = uicontrol('parent',  h.info_pnl_mesh,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              '120',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.12+0.315*(i-1) 0.34 0.2 0.165],...
                                     'HorizontalAlignment', 'center',...
                                     'visible',             'on');                                                                       
end                      
h.info_text_mesh_len = uicontrol('parent',  h.info_pnl_mesh,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'Mesh Length:',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0 0.166 0.32  0.165],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on');
cc = {'�_x:' '�_y:' '�_z:'};
for i = 1:3
    h.info_text_mesh_len_name(i) = uicontrol('parent',  h.info_pnl_mesh,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              cc{i},...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.05+0.315*(i-1) 0.01 0.07 0.165],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on'); 
    h.info_text_mesh_len_value(i) = uicontrol('parent',  h.info_pnl_mesh,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              '120',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.12+0.315*(i-1) 0.01 0.2 0.165],...
                                     'HorizontalAlignment', 'center',...
                                     'visible',             'on');                                                                       
end                     
%% Material Panel
h.info_text_material_type = uicontrol('parent',  h.info_pnl_material,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'Material Type:',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0 0.835 0.4 0.14],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on'); 
h.info_text_material_type_value = uicontrol('parent',  h.info_pnl_material,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'Isotropic',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.4 0.835 0.6 0.14],...
                                     'HorizontalAlignment', 'center',...
                                     'visible',             'on');
h.info_text_material_para = uicontrol('parent',  h.info_pnl_material,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'Material Parameters:',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0 0.67 0.4  0.165],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on');
cc  = {'�`:' '�g:' '�q:' '�^:'};                       
val = {'13' '1' '0' '0'};
for i = 1:4
    h.info_text_material_para_name(i) = uicontrol('parent',  h.info_pnl_material,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              cc{i},...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.01+0.245*(i-1) 0.505 0.07 0.165],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on'); 
    h.info_text_material_para_value(i) = uicontrol('parent',  h.info_pnl_material,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              val{i},...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.08+0.245*(i-1) 0.505 0.237 0.165],...
                                     'HorizontalAlignment', 'center',...
                                     'visible',             'on');                                                                       
end                                                      
h.info_text_material_file = uicontrol('parent',  h.info_pnl_material,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'Prototype:',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0 0.34 0.4 0.14],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on'); 
h.info_text_material_file_value = uicontrol('parent',  h.info_pnl_material,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'No002_Triclinic_P2I4',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.4 0.34 0.6 0.14],...
                                     'HorizontalAlignment', 'center',...
                                     'visible',             'on');
h.info_text_sphere_radius = uicontrol('parent',  h.info_pnl_material,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'Sphere Radius:',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0 0.17 0.4 0.14],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on'); 
h.info_text_sphere_radius_value = uicontrol('parent',  h.info_pnl_material,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              '0.12',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.4 0.17 0.6 0.14],...
                                     'HorizontalAlignment', 'center',...
                                     'visible',             'on');    
h.info_text_cylinder_radius = uicontrol('parent',  h.info_pnl_material,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              'Cylinder Radius:',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0 0.05 0.4 0.14],...
                                     'HorizontalAlignment', 'right',...
                                     'visible',             'on'); 
h.info_text_cylinder_radius_value = uicontrol('parent',  h.info_pnl_material,...
                                     'style',               'text',...
                                     'units',               'normalized',...
                                     'string',              '0.1',...
                                     'backgroundcolor',     bcolor,...
                                     'ForegroundColor',     [0 0 0],...
                                     'fontsize',            8,...
                                     'position',            [0.4 0.05 0.6 0.14],...
                                     'HorizontalAlignment', 'center',...
                                     'visible',             'on');                                    