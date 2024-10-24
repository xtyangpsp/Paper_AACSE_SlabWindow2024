function vardata = ncread_strict_dimension_order(ncFile, varName, dimNames, varargin)
%ncread_strict_dimension_order  Reads the data from the netcdf file in the specified dimension order
%
%   Function is as ncread, except that you can specify the dimension order of the output
%   This function is created tom solve the problem when reading netcdf
%   files generated by different engines like Matlab or snc.
%
%   The difference is that the variables can be written permutated.
%   Especially for 2D and 3D data this issue is important.
%
%   Syntax:
%       vardata = ncread_strict_dimension_order(ncFile, varName, dimNames)
%       vardata = ncread_strict_dimension_order(ncFile, varName, dimNames, start, count)
%       vardata = ncread_strict_dimension_order(ncFile, varName, dimNames, start, count, stride)
%
%   Input: 
%       ncFile    = url of the netcdf file to read
%       varName   = string with variable name to read
%       dimNames  = cell array with the names of the dimensions (strings)
%       $varargin = Extra parameters to pass to the ncread function like:
%                   start, count
%                   start, count, stride
%
%   Output:
%       vardata   = result of read data
%
%   Example
%       dimx = 'x_range';
%       dimy = 'y_range';
%       dimz = 'z_range';
%       V    = ncread_strict_dimension_order(ncFile,'volumedata', {dimx dimy dimz});
%
%
%   See also: ncread   

%% Copyright notice
%   --------------------------------------------------------------------
%   Copyright (C) 2012 Van Oord
%       Thijs Damsma
%
%       tda@vanoord.com
%
%       Watermanweg 64
%       3067 GG
%       Rotterdam
%       Netherlands
%
%   This library is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
%
%   This library is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You should have received a copy of the GNU General Public License
%   along with this library.  If not, see <http://www.gnu.org/licenses/>.
%   --------------------------------------------------------------------

% This tool is part of <a href="http://www.OpenEarth.eu">OpenEarthTools</a>.
% OpenEarthTools is an online collaboration to share and manage data and
% programming tools in an open source, version controlled environment.
% Sign up to recieve regular updates of this function, and to contribute
% your own tools.

%% Version <http://svnbook.red-bean.com/en/1.5/svn.advanced.props.special.keywords.html>
% Created: 09 Oct 2012
% Created with Matlab version: 8.0.0.783 (R2012b)

% $Id$
% $Date$
% $Author$
% $Revision$
% $HeadURL$
% $Keywords: $
assert(iscellstr(dimNames),'dimNames must be a cell string')

info = ncinfo(ncFile, varName);
assert(numel(dimNames) == numel(info.Dimensions),'Wrong number of dimensions')

[in_a,order] = ismember({info.Dimensions.Name},dimNames);
assert(all(in_a),'Wrong dimension names, should be these');
% are all
if issorted(order)
    vardata = ncread(ncFile, varName, varargin{:});
else
    % Permute start, count, stride
    varargin = cellfun(@(arg) arg(order),varargin,'UniformOutput',false);
    varargin
    % read data and permute outcome
    vardata = permute(ncread(ncFile, varName, varargin{:}),order);
end