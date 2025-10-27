function mustBeCell(A)
%MUSTBECELL Validate that value is cell
%   MUSTBECELL(A) throws an error if A is not a cell array.
%   MATLAB call iscell to determine if a value is a cell.
%
%   See also: ISCELL.

if ~isempty(A) & ~iscell(A)
    throwAsCaller(MException('sayHello:inputError','Input must be cell or empty.'));
end
