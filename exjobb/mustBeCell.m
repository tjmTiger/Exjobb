function mustBeCell(A)
%MUSTBENUMERIC Validate that value is cell
%   MUSTBENUMERIC(A) throws an error if A contains nonnumeric values.
%   MATLAB call iscell to determine if a value is a cell.
%
%   See also: ISCELL.

if ~isempty(A) & ~iscell(A)
    throwAsCaller(MException('sayHello:inputError','Input must be cell or empty.'));
end
