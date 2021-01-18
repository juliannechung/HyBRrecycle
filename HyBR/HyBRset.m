function options = HyBRset(varargin)
% function options = HyBRset(varargin)
%
%   OPTIONS = HyBRset(varargin)
%
%   Create/alter options structure for HyBR code.
%   OPTIONS = HyBRset('PARAM1',VALUE1,'PARAM2',VALUE2,...) 
%     creates an options structure in which the named parameters have
%     the specified values.  Any unspecified parameters are set to [] 
%     (parameters with value [] indicate to use the default value for that 
%     parameter when passed to the HyBR function). It is sufficient to type
%     only the leading characters that uniquely identify the parameter.  
%     Case is ignored for parameter names.
%     NOTE: For values that are strings, the complete string is required.
%
%   OPTIONS = HyBRset(OLDOPTS,'PARAM1',VALUE1,...) creates a copy of
%     OLDOPTS with the named parameters altered with the specified values.
%
%   OPTIONS = HyBRset(OLDOPTS,NEWOPTS) combines an existing options structure
%     OLDOPTS with a new options structure NEWOPTS.  Any parameters in NEWOPTS
%     with non-empty values overwrite the corresponding old parameters in
%     OLDOPTS.
%
%   HyBRset with no input arguments and no output arguments displays all
%     parameter names and their possible values, with defaults shown in {}
%
%   OPTIONS = HyBRset(with no input arguments) creates an options structure
%     where all the fields are set to [].
%
%   OPTIONS = HyBRset('HyBR') creates an options structure with all
%     the parameter names and default values relevant to 'HyBR'. That is,
%           HyBRset('HyBR')
%   or
%           HyBRset(@HyBR)
%   returns an options structure containing all the parameter names and
%   default values relevant to the function 'HyBR'.
%
% HyBRset PARAMETERS for MATLAB: ( default parameter in {} )
%   InSolv - Solver for the inner problem: [ none | TSVD | {Tikhonov} ]
%   RegPar - Either (a) a value for the regularization parameter
%                           [ real non-negative scalar ]
%                   (b) a method for choosing a reg. parameter
%                           [ DP | GCV | {WGCV} ] 
%                       where  (i) DP - discrepancy principle
%                              (ii) 'GCV' - standard GCV
%                              (iii) 'WGCV' - weighted GCV
%                   (c) finds the optimal reg. parameter: 
%                           [ optimal ]     (requires x_true)
%   nLevel - if RegPar is 'DP' or 'UPRE', then nLevel represents the noise level and must be
%                       [non-negative scalar | {est}]
%   Omega - If RegPar = 'WGCV', then omega must be either:
%                           [ non-negative scalar | {adapt} ] 
%            where  (a) non-negative scalar - a value for omega
%                   (b) adapt - uses the adaptive method
%   Iter - Maximum number of Lanczos iterations 
%                   [positive integer | {[]} ]
%   Reorth - Reorthogonalize Lanczos subspaces: [ on | {off} ]
%   x_true - True solution : [ array | {off} ]
%                Returns error with respect to x_true at each iteration
%                and is used to compute 'optimal' regularization parameters
%   BegReg - Begin regularization after this iteration: 
%                [ positive integer | {2} ]
%  FlatTol - Tolerance for detecting flatness in the GCV curve as a
%               stopping criteria
%               [ non-negative scalar | {10^-6}]
%   MinTol - Window of iterations for detecting a minimum of the GCV curve
%               as a stopping criteria
%               [ positive integer | {3}]
%   ResTol - Residual tolerance for stopping the LBD iterations
%                   [[non-negative scalar, non-negative scalar]  | {[10^-6, 10^-6]}]
%
%   Examples:
%     To create OPTIONS with the default options for HyBR
%       OPTIONS = HyBRset('HyBR');
%     To create an OPTIONS structure with RegPar = 'WGCV' and Omega = .5
%       OPTIONS = HyBRset('RegPar','WGCV', 'Omega',.5);
%     To change the maximum iterations to 150 in OPTIONS
%       OPTIONS = HyBRset(OPTIONS,'Iter',150);
%
%   See also HyBR.
%
%   J.Chung 3/2014
%

% Print out possible values of properties. 
if (nargin == 0) && (nargout == 0)
    fprintf('            InSolv: [ none | TSVD | {Tikhonov} ]\n');
    fprintf('            RegPar: [ non-negative scalar | DP | GCV | {WGCV} | upre| optimal|NWGCV ]\n');
    fprintf('            nLevel: [ non-negative scalar | {est} ]\n');
    fprintf('             Omega: [ non-negative scalar | {adapt} ]\n');
    fprintf('              Iter: [ positive integer  | {[]} ]\n');
    fprintf('            Reorth: [ on | {off} ]\n');
    fprintf('            x_true: [ array | {off} ]\n');
    fprintf('            BegReg: [ positive integer | {2} ]\n');
    fprintf('           FlatTol: [  non-negative scalar | {10^-6}  ]\n');
    fprintf('            MinTol: [  positive integer | {3}  ]\n');
    fprintf('            ResTol: [  non-negative scalar | {[10^-6 10^-6]}  ]\n');
    return;
end

% Create a struct of all the fields
allfields = {'InSolv'; 'RegPar';'nLevel';'Omega';'Iter';'Reorth'; ...
    'x_true';'BegReg'; 'FlatTol'; 'MinTol'; 'ResTol'};
  
% create cell array
structinput = cell(2,length(allfields));
% fields go in first row
structinput(1,:) = allfields';
% []'s go in second row
structinput(2,:) = {[]};
% turn it into correctly ordered comma separated list and call struct
options = struct(structinput{:});

numberargs = nargin; % we might change this value, so assign it


% If we pass in a function name then return the defaults.
if (numberargs==1) && (ischar(varargin{1}) || isa(varargin{1},'function_handle') )
    if ischar(varargin{1})
        funcname = lower(varargin{1});
        if ~exist(funcname)
            error('Undefined function.  Please use HyBRset(@HyBR).')
        end
    elseif isa(varargin{1},'function_handle')
        funcname = func2str(varargin{1});
    end
    try 
      optionsfcn = feval(varargin{1},'defaults');
    catch
        error('HyBRset ONLY works with HyBR.  Please use HyBRset(@HyBR).')
    end
    % The defaults from the HyBR functions don't include all the fields
    % typically, so run the rest of HyBRset as if called with
    % HyBRset(options,optionsfcn) to get all the fields.
    varargin{1} = options;
    varargin{2} = optionsfcn;
    numberargs = 2;
end

Names = allfields;
m = size(Names,1);
names = lower(Names);

i = 1;
while i <= numberargs
    arg = varargin{i};
    if ischar(arg)                         % arg is an option name
        break;
    end
    if ~isempty(arg)                      % [] is a valid options argument
        if ~isa(arg,'struct')
            error(['Expected argument %d to be a string parameter name ' ...
                'or an options structure \n created with HyBRset.'], i);
        end
        for j = 1:m
            if any(strcmp(fieldnames(arg),Names{j,:}))
                val = arg.(Names{j,:});
            else
                val = [];
            end
            if ~isempty(val)
                if ischar(val)
                    val = lower(deblank(val));
                end
                checkfield(Names{j,:},val)
                options.(Names{j,:}) = val;
            end
        end
    end
    i = i + 1;
end

% A finite state machine to parse name-value pairs.
if rem(numberargs-i+1,2) ~= 0
    error('Arguments must occur in name-value pairs.');
end

expectval = 0;                          % start expecting a name, not a value
while i <= numberargs
    arg = varargin{i};

    if ~expectval
        if ~ischar(arg)
            error('Expected argument %d to be a string parameter name.', i);
        end

        lowArg = lower(arg);
        j = strmatch(lowArg,names);
        if isempty(j)                       % if no matches
            error( sprintf('Invalid parameter name ''%s'' ', arg));
        elseif length(j) > 1                % if more than one match
            % Check for any exact matches (in case any names are subsets of others)
            k = strmatch(lowArg,names,'exact');
            if length(k) == 1
                j = k;
            else
                error(sprintf('Ambiguous parameter name ''%s'' ', arg));
            end
        end
        expectval = 1;                      % we expect a value next

    else
        if ischar(arg)
            arg = lower(deblank(arg));
        end
        checkfield(Names{j,:},arg);
        options.(Names{j,:}) = arg;
        expectval = 0;
    end
    i = i + 1;
end

if expectval
    error( 'Expected value for parameter ''%s''.', arg);
end


%----SUBFUNCTION---------------------------------------------
function checkfield(field,value)
%CHECKFIELD Check validity of structure field contents.
%   CHECKFIELD('field',V) checks the contents of the specified
%   value V to be valid for the field 'field'. 
%

% empty matrix is always valid
if isempty(value)
    return
end

% See if 'field' is a valid field.
validfield = true;
switch field
    case {'InSolv'} % none, tsvd, tikhonov
        [validvalue, errmsg] = InSolvetype(field,value);
    case {'RegPar'} % real non-negative scalar, gcv, wgcv, optimal
        [validvalue, errmsg] = RegPartype(field,value);
    case {'nLevel'} % real non-negative scalar, adapt
        [validvalue, errmsg] = nLeveltype(field,value);  
    case {'Omega'} % real non-negative scalar, adapt
        [validvalue, errmsg] = Omegatype(field,value);
    case {'Iter'} % real positive integer
        [validvalue, errmsg] = PosInteger(field,value);
    case {'Reorth'} % off, on
        [validvalue, errmsg] = onOffType(field,value);
    case {'x_true'} % numeric array, off
        [validvalue, errmsg] = x_truetype(field,value);
    case {'BegReg'}% real positive integer
        [validvalue, errmsg] = PosInteger(field,value);
    case {'FlatTol'}% real non-negative scalar
        [validvalue, errmsg] = nonNegscalar(field,value);
    case {'MinTol'}% real positive integer
        [validvalue, errmsg] = PosInteger(field,value);
    case {'ResTol'}% real non-negative vector of length 2
        [validvalue, errmsg] = nonNeg2vector(field,value);
  otherwise
        validfield = false;  
        validvalue = false;
        errmsg = sprintf('Unrecognized parameter name ''%s''.', field);
end

if validvalue 
    return;
else
  error(errmsg)
end

%------------------------------------------------------------------------

function [valid, errmsg] = InSolvetype(field,value)
% One of these strings: tsvd, tikhonov, none
valid =  ischar(value) && any(strcmp(value,{'tsvd','tikhonov','none'}));
if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''tsvd'', ''tikhonov'' or ''none''.',field);
else
  errmsg = '';
end

%--------------------------------------------------------------------------

function [valid, errmsg] = RegPartype(field,value)
% One of these: real nonnegative scalar, DP, GCV, WGCV, optimal
valid =  (isreal(value) && isscalar(value) && (value >= 0)) | (ischar(value) && any(strcmp(value,{'dp', 'gcv', 'wgcv', 'upre','optimal','nwgcv'})));

if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a non-negative scalar or ''DP'' or ''GCV'' or ''WGCV'' or ''upre'' or''optimal'' or ''NWGCV''.',field);
else
  errmsg = '';
end

%--------------------------------------------------------------------------

function [valid, errmsg] = nLeveltype(field,value)
% One of these: real non-negative scalar, adapt
valid =  (isreal(value) && isscalar(value) && (value >= 0) )| (ischar(value) && any(strcmp(value,{'est'})));

if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a non-negative scalar or ''est''.',field);
else
  errmsg = '';
end

%--------------------------------------------------------------------------

function [valid, errmsg] = Omegatype(field,value)
% One of these: real non-negative scalar, adapt
valid =  (isreal(value) && isscalar(value) && (value >= 0) )| (ischar(value) && any(strcmp(value,{'adapt'})));

if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a non-negative scalar or ''adapt''.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = PosInteger(field,value)
% Any positive real integer
valid =  isreal(value) && isscalar(value) && (value > 0) && value == floor(value) ;

if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a positive real integer.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = onOffType(field,value)
% One of these strings: on, off
valid =  ischar(value) && any(strcmp(value,{'on';'off'}));
if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be ''off'' or ''on''.',field);
else
  errmsg = '';
end

%-----------------------------------------------------------------------

function [valid, errmsg] = x_truetype(field,value)
% Either a numeric array, or off
valid =  isnumeric(value) | (ischar(value) && any(strcmp(value,{'off'})));
if ~valid
    errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a numeric array or ''off''.',field);
else
  errmsg = '';
end


%------------------------------------------------------------------------

function [valid, errmsg] = nonNegscalar(field,value)
% Any real non-negative scalar
valid =  isreal(value) && isscalar(value) && (value >= 0);

if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative scalar.',field);
else
  errmsg = '';
end

%------------------------------------------------------------------------

function [valid, errmsg] = nonNeg2vector(field,value)
% Any real non-negative scalar
valid =  isreal(value) && isvector(value) && sum(value>=0)==length(value) && (length(value)==2);

if ~valid
 errmsg = sprintf('Invalid value for OPTIONS parameter %s: must be a real non-negative vector of length 2.',field);
else
  errmsg = '';
end
