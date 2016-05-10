% Matrix Expression Class for Automatic Differentiation
% Emanuele Ruffaldi 2016
%
% Provides Reverse-Mode Automatic Differentiation of Matrices
classdef matexp < handle
    
    properties
        aname
        avalue
        aadjoint
        aop
        aoperands
    end
    
    methods
        % Possible constructors:
        % matexp(value)
        % matexp(name,value)
        % matexp([],op,operands)
        function this = matexp(varargin)
            if nargin == 2
                this.aname = varargin{1};
                this.avalue = varargin{2};
                this.aop= '';
            elseif nargin == 1
                this.avalue = varargin{1};
                this.aop= '';
            elseif nargin == 3
                % skip first
                this.aop = varargin{2};
                this.aoperands = varargin{3};
            else
                error('Invalid matexp constructor');
            end
        end

        % computes autodiff in recursive mode 
        % ASSUMES an update of the value
        function autodiff(this,isroot)
            if nargin == 1
                % the root performs the reset of the adjoints
                resetadjoint(this,0);
                this.aadjoint = eye(length(this.avalue));
            end
            
            ops = this.aoperands;
            A = this.aadjoint;
            V = this.avalue;
            
            % by operation, increment adjoint of children
            switch(this.aop)
                case '+'
                    ops{1}.aadjoint = ops{1}.aadjoint + A;
                    ops{2}.aadjoint = ops{2}.aadjoint + A;
                case '-'
                    ops{1}.aadjoint = ops{1}.aadjoint + A;
                    ops{2}.aadjoint = ops{2}.aadjoint - A;
                case '*'
                    ops{1}.aadjoint = ops{1}.aadjoint + ops{2}.avalue*A;
                    ops{2}.aadjoint = ops{2}.aadjoint + A*ops{1}.avalue;
                case 'logdet'
                    ops{1}.aadjoint = ops{1}.aadjoint + inv(v)'; % X^(-T)
                case 'det'
                    assert('Not implemented autodiff of det');
                case 't' 
                    ops{1}.aadjoint = A; % propagate TODO
                case 'vec'
                case 'inv'
                    ops{1}.aadjoint = ops{1}.aadjoint - V*A*V;
                case 'tr'
                    ops{1}.aadjoint = A'; % propagate transpose
                case ''  % nothing
                    return
                otherwise
                    error(['Unimplemented ' this.aop]);
            end
            this.aoperands = ops;
            
            % then continue the descent
            for I=1:length(this.aoperands)
                autodiff(this.aoperands{I},0);
            end
        end
        
        
        % collects the variables present in the expression tree
        function r = collectvars(this)
            r = {};
            if  ~isempty(this.aname)
                r = {this};
            else
                for I=1:length(this.aoperands)
                    w = collectvars(this.aoperands{I});
                    r = [r; w];
                end
            end            
        end

        % resets he adjoint before autodiff
        function resetadjoint(this,x)
            for I=1:length(this.aoperands)
                resetadjoint(this.aoperands{I},x);
            end
            this.aadjoint = x;
        end

        % updates the value
        function update(this)
            for I=1:length(this.aoperands)
                update(this.aoperands{I});
            end
            switch(this.aop)
                case '+'
                    this.avalue = this.aoperands{1}.avalue+this.aoperands{2}.avalue;
                case '-'
                    this.avalue = this.aoperands{1}.avalue-this.aoperands{2}.avalue;
                case '*'
                    this.avalue = this.aoperands{1}.avalue*this.aoperands{2}.avalue;
                case 'inv'
                    this.avalue = inv(this.aoperands{1}.avalue);
                case 'tr'
                    this.avalue = (this.aoperands{1}.avalue)';
                case 't'
                    this.avalue = trace(this.aoperands{1}.avalue);
                case 'logdet'
                    this.avalue = log(det(this.aoperands{1}.avalue));
                case 'det'
                    this.avalue = det(this.aoperands{1}.avalue);
                case 'vec'
                    c = this.aoperands{1}.avalue;
                    this.avalue = c(:);                    
                case '' % leaf                   
                otherwise
                    error('unsupported operator');
            end
        end
                
        function r = plus(a,b)
            if ~isa(b,'matexp')
                b = matexp(b);
            end
            if ~isa(a,'matexp')
                a = matexp(a);
            end
            r = matexp([],'+',{a,b});
        end        
        
        function r = minus(a,b)
            if ~isa(b,'matexp')
                b = matexp(b);
            end
            if ~isa(a,'matexp')
                a = matexp(a);
            end
            r = matexp([],'-',{a,b});
        end        
        
        function r = mtimes(a,b)
            if ~isa(b,'matexp')
                b = matexp(b);
            end
            if ~isa(a,'matexp')
                a = matexp(a);
            end
            r = matexp([],'*',{a,b});
        end        
        
        function r = ctranspose(a)
            r = matexp([],'tr',{a});
        end        
        
        function r = trace(a)
            r = matexp([],'t',{a});
        end        
        
        function r = inv(this)
            r = matexp([],'inv',{this});
        end

        function r = det(this)
            r = matexp([],'det',{this});
        end
                
        function r = log(this)
            assert(strcmp(this.aop,'det'),'Only log det supported');
            r = matexp([],'logdet',{this.aoperands{1}});
        end
        
        % make column vector
        function r = vec(this)
            r = matexp([],'vec',{this});
        end
        
        % size of the value
        function r = size(this)
            r = size(this.avalue);
        end
        
        % returns the value
        function r = value(this)
            r = this.avalue;
        end
        
        % name for variables
        function r = name(this)
            r = this.aname;
        end
        
        % returns adjoint
        function r = adjoint(this)
            r = this.aadjoint;
        end
        
        % sets the value for the constant or variable ones
        function set(this,value)
            assert(isempty(this.aop));
            this.avalue = value;
        end         
    end
    
end

