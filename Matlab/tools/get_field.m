function [field1,field2,field3,field4,field5,field6,field7,field8,field9,field10] = get_field(file, varnam_cell)

% GET_AVG_FIELD_3d  Reads fields from a set of netCDF files and averages them.
%  
%   FIELD = GET_AVG_FIELD(FILES, VARNAM) reads the fields with
%   variable name VARNAM (can be a cell array of variable names)
%   from the netCDF files given in the cell FILES. The average 
%   of these fields is returned as FIELD1-FIELD10. Output is only 
%   given up to the nuber of variables specified in VARNAM. 
%   Maximum number of variables is 10

%   Output format is in the format the variables are stored in netcdf

field1 = 0;
field2 = 0;
field3 = 0;
field4 = 0;
field5 = 0;
field6 = 0;
field7 = 0;
field8 = 0;
field9 = 0;
field10 = 0;

if iscell(varnam_cell)
for i = 1:length(varnam_cell) % loop over variables in VARNAM
    varnam = varnam_cell{i};
  
  %find alias for pre-fram simulations.
 % varname = alias(files, varnam);
  varname = varnam;
  
  % read fields from files
  ncid     = netcdf.open(char(file), 'NC_NOWRITE');
  field    = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,char(varname))));

    netcdf.close(ncid);  
  
switch i
    case 1
        field1 = field;
    case 2
        field2 = field;
    case 3
        field3 = field;
    case 4
        field4 = field;
    case 5
        field5 = field;
    case 6
        field6 = field;
    case 7
        field7 = field;
    case 8
        field8 = field;
    case 9
        field9 = field;
    case 10
        field10 = field;
    otherwise
        disp('get_avg_field gets up to 10 variables, all subsequent ignored')
end
  
end
else
  varnam = varnam_cell;
  
  %find alias for pre-fram simulations.
 % varname = alias(files, varnam);
  varname = varnam;
  
  % read fields from files

    ncid     = netcdf.open(char(file), 'NC_NOWRITE');
    field    = double(netcdf.getVar(ncid,netcdf.inqVarID(ncid,char(varname))));  
    netcdf.close(ncid);  

  field1 = field;
end
      