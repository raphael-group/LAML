class Param():
    def __init__(self,names,values,lower_bounds,upper_bounds):
        self.names = names
        self.values = values
        self.lower_bounds = lower_bounds
        self.upper_bounds = upper_bounds
        self.name2bound = {n:(l,u) for (n,l,u) in zip(names,lower_bounds,upper_bounds)}

    def get_names(self):
        return self.names

    def get_bound(self,pname):
        return self.name2bound[pname]    

    def get_name2value_dict(self):
        return {n:v for n,v in zip(self.names,self.values)}

    def get_value(self,pname):
        for i,(n,v) in enumerate(zip(self.names,self.values)):
            if n == pname:
                return v
        return None # only get here if pname is not found in self.names    

    def show_values(self):
        out_str = ""
        for n,v in zip(self.names,self.values):
            out_str += str(n) + "=" + str(v) + " "
        return out_str    

    def set_value(self,pname,value):
        for i,n in enumerate(self.names):
            if n == pname:
                self.values[i] = value
                return True 
        return False # only get here if pname is not found in self.names    

