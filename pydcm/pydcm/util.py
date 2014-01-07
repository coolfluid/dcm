import coolfluid as cf
import pickle

KEY=0; VAL=1;
fparser_mapping = {
      '[0]':'__0', '[1]':'__1', '[2]':'__2', '[3]':'__3', '[4]':'__4',
      '[x]':'__0', '[y]':'__1', '[z]':'__2' }
eval_mapping = {
      '[x]':'[0]', '[y]':'[1]', '[z]':'[2]' }
      
def to_eval_str(arg):
    if isinstance(arg,str) or isinstance(arg,unicode):
        eval_str = arg.replace('^','**')
        for k, v in eval_mapping.iteritems():
            eval_str = eval_str.replace(k, v)
        return eval_str
    else:
        eval_str = 'float('+str(arg)+')'
        return eval_str
        
def to_fparser_str(arg):
    if isinstance(arg,str) or isinstance(arg,unicode):
        fparser_str = arg.replace('**','^')
        for k, v in fparser_mapping.iteritems():
            fparser_str = fparser_str.replace(k, v)
        return fparser_str
    else:
        fparser_str = str(arg)
        return fparser_str
    
class Definitions(object):
    """
    Definitions is a helper class that allows
    """
    
    def __init__(self):
        self.expressions = []
        self.position = {}

    def set(self,**kwargs):
        for (key,value) in kwargs.iteritems():
            if key in self.position:
                self.expressions[self.position[key]][VAL] = value
            else:
                self.position[key] = len(self.expressions)
                self.expressions.append( [key,value] )
                
    def has(self,var):
        return (var in self.position)

    def fparser(self,expr):
        expressions = []
        for (key, val) in self.expressions:
            if not isinstance(val, list):
                expressions.append(str(key)+':='+to_fparser_str(val))
            else:
                for n,v in enumerate(val):
                    expressions.append(str(key)+'__'+str(n)+':='+to_fparser_str(v))
        expressions_str = '; '.join(expressions)+'; '+to_fparser_str(expr)
        return expressions_str

    def eval(self,var):
        from math import sqrt, log, exp, sin, cos, tan
        pos = len(self.expressions)
        if self.has(var):
             pos = self.position[var]
        for expr in self.expressions[:pos]:
            if isinstance(expr[VAL],list):
                exec(expr[KEY]+'=[]')
                for v in expr[VAL]:
                    exec(expr[KEY]+'.append('+to_eval_str(v)+')')
            else:
                exec(expr[KEY]+' = '+to_eval_str(expr[VAL]))

        if self.has(var):
            expr = self.expressions[pos]
            if isinstance(expr[VAL],list):
                exec(expr[KEY]+'=[]')
                for v in expr[VAL]:
                    exec(expr[KEY]+'.append('+to_eval_str(v)+')')
                return eval(expr[KEY])
            else:
                return eval(to_eval_str(expr[VAL]))
        else:
            return eval(var)
            
def log(*args):
    if cf.Core.rank() == 0:
        print(" ".join([str(arg) for arg in args]))
    
def write(obj,filename):
    if cf.Core.rank() == 0:
        pickle.dump( obj, open(filename, 'w') )        
    
def read(filename):
    return pickle.load( open( filename, 'r' ) )
    
def symlink(src,dest):
    if cf.Core.rank() == 0:
        import os, errno
        try: 
            os.symlink(src, dest) 
        except OSError, e: 
            if e.errno == errno.EEXIST: 
                os.remove(dest) 
                os.symlink(src, dest)
            
def mkdir(dir):
    if cf.Core.rank() == 0:
        import os, errno 
        try: 
            os.mkdir(dir) 
        except OSError, e: 
            if e.errno == errno.EEXIST: 
                return
            raise OSError,str(e)