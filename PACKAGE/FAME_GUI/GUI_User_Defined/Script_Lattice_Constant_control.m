h.text_constant_a = uicontrol( 'parent',          h.frame_constant,...
                               'style',           'text',...
                               'units',           'normalized',...
                               'string',          'a = ',...
                               'backgroundcolor', bcolor,...
                               'ForegroundColor', [0 0 0],...
                               'fontsize',        12,...
                               'position',        [.08 .84 .1 0.14]);
h.edit_constant_a = uicontrol( 'parent',          h.frame_constant,...
                               'style',           'edit',...
                               'units',           'normalized',...
                               'ForegroundColor', [0 0 0],...
                               'string',          '1',...
                               'fontsize',        12,...
                               'position',        [.19 .84 .25 0.14]  );
h.text_constant_b = uicontrol( 'parent',          h.frame_constant,...
                               'style',           'text',...
                               'units',           'normalized',...
                               'string',          'b = ',...
                               'backgroundcolor', bcolor,...
                               'ForegroundColor', [0 0 0],...
                               'fontsize',        12,...
                               'position',        [.08 0.68 .1 0.14]);
h.edit_constant_b = uicontrol( 'parent',          h.frame_constant,...
                               'style',           'edit',...
                               'units',           'normalized',...
                               'ForegroundColor', [0 0 0],...
                               'string',          '1',...
                               'fontsize',        12,...
                               'position',        [.19 .68 .25 0.14],...
                               'enable', 'off');
h.text_constant_c = uicontrol( 'parent',          h.frame_constant,...
                               'style',           'text',...
                               'units',           'normalized',...
                               'string',          'c = ',...
                               'backgroundcolor', bcolor,...
                               'ForegroundColor', [0 0 0],...
                               'fontsize',        12,...
                               'position',        [.08 0.52 .1 0.14]);
h.edit_constant_c = uicontrol( 'parent',          h.frame_constant,...
                               'style',           'edit',...
                               'units',           'normalized',...
                               'ForegroundColor', [0 0 0],...
                               'string',          '1',...
                               'fontsize',        12,...
                               'position',        [.19 .52 .25 0.14],...
                               'enable', 'off'); 
h.text_constant_alpha = uicontrol( 'parent',          h.frame_constant,...
                               'style',           'text',...
                               'units',           'normalized',...
                               'string',          '�\ = ',...
                               'backgroundcolor', bcolor,...
                               'ForegroundColor', [0 0 0],...
                               'fontsize',        12,...
                               'position',        [.5 .84 .1 0.14]);
h.edit_constant_alpha = uicontrol( 'parent',          h.frame_constant,...
                               'style',           'edit',...
                               'units',           'normalized',...
                               'ForegroundColor', [0 0 0],...
                               'string',          '1',...
                               'fontsize',        12,...
                               'position',        [.61 .84 .25 0.14],...
                               'enable', 'off');
h.text_constant_beta = uicontrol( 'parent',          h.frame_constant,...
                               'style',           'text',...
                               'units',           'normalized',...
                               'string',          '�] = ',...
                               'backgroundcolor', bcolor,...
                               'ForegroundColor', [0 0 0],...
                               'fontsize',        12,...
                               'position',        [.5 0.68 .1 0.14]);
h.edit_constant_beta = uicontrol( 'parent',          h.frame_constant,...
                               'style',           'edit',...
                               'units',           'normalized',...
                               'ForegroundColor', [0 0 0],...
                               'string',          '1',...
                               'fontsize',        12,...
                               'position',        [.61 .68 .25 0.14],...
                               'enable', 'off');
h.text_constant_gamma = uicontrol( 'parent',          h.frame_constant,...
                               'style',           'text',...
                               'units',           'normalized',...
                               'string',          '�^ = ',...
                               'backgroundcolor', bcolor,...
                               'ForegroundColor', [0 0 0],...
                               'fontsize',        12,...
                               'position',        [.5 0.52 .1 0.14]);
h.edit_constant_gamma = uicontrol( 'parent',          h.frame_constant,...
                               'style',           'edit',...
                               'units',           'normalized',...
                               'ForegroundColor', [0 0 0],...
                               'string',          '1',...
                               'fontsize',        12,...
                               'position',        [.61 .52 .25 0.14],...
                               'enable', 'off');
h.frame_text_information = uipanel( 'parent',          h.frame_constant,...
                                          'units',           'normalized',...
                                          'backgroundcolor', bcolor,...
                                          'ForegroundColor', [0 0 0],...
                                          'shadowcolor',     [0 0 0],...
                                          'title',           'Notice',...
                                          'position',        [.05 .01 .9 0.5],... 
                                          'fontsize',        12,...
                                          'visible',         'on');                  

h.text_information    = uicontrol( 'parent',          h.frame_text_information,...
                               'style',           'text',...
                               'units',           'normalized',...
                               'string',          'lattice constant must satisfying a = b = c, �\=�]=�^=�k/2',...
                               'HorizontalAlignment', 'left',...
                               'backgroundcolor', bcolor,...
                               'ForegroundColor', [0 0 0],...
                               'fontsize',        12,...
                               'position',        [.01 .01 .98 0.98]);