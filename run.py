"""
MIT License

Copyright (c) 2017 Joris Gillis, KU Leuven

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""

import sys
import os
#import ipdb
from pycparser import c_parser, c_ast, parse_file, c_generator, CParser
from collections import defaultdict 

import os

basepath = os.path.dirname(os.path.realpath(__file__))

model_file_name = "Drivetrain_Rigid_Model"

if len(sys.argv)>1:
  model_file_name = sys.argv[1]

class Sparsity:
  def __init__(self):
    pass
    
  def pack(self):
    pass
    
  def __repr__(self):
    return str((self.rows, self.cols))
    
  def colind(self):
    if hasattr(self,"mJc"):
      colind = [0]*(self.cols+1)
      for k,v in self.mJc:
        colind[k] = v
      return colind
    else:
      return None

  def row(self):
    if hasattr(self,"mIr"):
      row = [0]*len(self.mIr)
      for k,v in self.mIr:
        row[k] = v
      return row
    else:
      return None
      
matrices = defaultdict(Sparsity)

codes = {}
codes_nodes = {}

code_names = ["m","a","b","y","dxy","f","tduf","tdxf","r","mode","dp_r","del_v","del_v0","del_tmax","del_t"]



#   m*dot(x) = a*x + b*u + f(x,u)
#      y(x)
# dxy:  dy/dx
# dxf:  df/dx
# tdxf: ? -- has no method to specify nonzero values, even though sparsity has nonzeros
# tduf: ? -- has no method to specify nonzero values, even though sparsity has nonzeros
# z : ?
# r : ?

# Number of dynamic variable constraints: !0 -> higher index DAE
# Can see eliminated vars through statistics viewer

#Simulink.BlockDiagram.saveActiveConfigSet('Drivetrain_Rigid_Model', 'my_config_set.m')

#https://nl.mathworks.com/help/symbolic/reduceredundancies.html

"""

  }
  ContStates {
    ContStateDefaults {
      RecordType	      ContState
      Width		      1
      StartIndex	      -1
      DataLoggingOn	      1
      InitialValue	      [0.0]
      StorageClass	      "Auto"
      HasObject		      0
      MemoryMapIdx	      [-1, -1, -1]
      SigSrc		      [-1, -1, -1]
      DerivativeMemoryMapIdx  [-1, -1, -1]
    }
    NumContStates	    1
    ContState {
      VarGroupIdx	      [2, 3, 4, -1, 0]
      Identifier	      "Drivetrain_Rigid_ModelFlywheel_"
      Width		      2
      StartIndex	      0
      GrSrc		      [7, 4]
      Partitions {
	NumPartitions		2
	Partition {
	  Name			  "Drivetrain_Rigid_Model.Flywheel_Inertia.w"
	  PathAlias		  "Drivetrain_Rigid_Model/Flywheel Inertia"
	  Width			  1
	  CStateIndices		  [0]
	}
	Partition {
	  Name			  "Drivetrain_Rigid_Model.MotionSensor.Ideal_Rotational_Motion_Sensor.phi"
	  PathAlias		  "Drivetrain_Rigid_Model/MotionSensor/Ideal Rotational Motion Sensor"
	  Width			  1
	  CStateIndices		  [1]
	}
      }
      
"""

class PreTransform(c_ast.NodeVisitor):

  def visit_Compound(self, n):
      block_items = n.block_items or []
      
      for i,b in enumerate(block_items):
        if isinstance(b,c_ast.UnaryOp) and b.op=="p--":
          block_items[i] = c_ast.Assignment("-=", b.expr, c_ast.Constant("int","1"))
        elif isinstance(b,c_ast.UnaryOp) and b.op=="p++":
          block_items[i] = c_ast.Assignment("+=", b.expr, c_ast.Constant("int","1"))
        else:
          self.visit(b)

  
  
def pre_transform(a):
  v = PreTransform()
  v.visit(a)
  return a

generator = c_generator.CGenerator()
def to_c(node):
  return generator.visit(node)

def from_constant(n):
     if n.value.endswith("UL"):
       return n.value[:-2]
     elif n.value.endswith("ULL"):
       return n.value[:-3]
     else:
       return n.value


class CommonExpressionGenerator(c_generator.CGenerator):
   def visit_Constant(self, n):
     return from_constant(n)
    
class MatlabExpressionGenerator(CommonExpressionGenerator):
  def visit_ArrayRef(self, n):
    arrref = self._parenthesize_unless_simple(n.name)
    return arrref + '(' + self.visit(n.subscript) + '+1)'
    
  def visit_StructRef(self, n):
      sref = self.visit(n.name)
      return sref + "." + self.visit(n.field)

  def visit_UnaryOp(self, n):
    operand = self._parenthesize_unless_simple(n.expr)
    if n.op=="&":
      return operand
    if n.op == '!':
      return '~' + operand
    else:
      return CommonExpressionGenerator.visit_UnaryOp(self, n)

  def visit_BinaryOp(self, n):
    m = {"!=": "~=", "&&": "&", "||": "|"}
    if n.op in m.keys():
        lval_str = self._parenthesize_if(n.left,
                            lambda d: not self._is_simple_node(d))
        rval_str = self._parenthesize_if(n.right,
                            lambda d: not self._is_simple_node(d))
        return '%s %s %s' % (lval_str, m[n.op], rval_str)
    else:
      return CommonExpressionGenerator.visit_BinaryOp(self, n)

  def visit_ID(self, n):
      return self.to_matlab_name(n.name)

  def visit_Compound(self, n):
      s = self._make_indent() + '\n'
      self.indent_level += 2
      if n.block_items:
          s += ''.join(self._generate_stmt(stmt) for stmt in n.block_items)
      self.indent_level -= 2
      s += self._make_indent() + '\n'
      return s
      
  def to_matlab_name(self, name):
      if name.startswith('_'):
        return "d" + name
      else:
        return name

  def visit_Assignment(self, node):
      rval_str = self._parenthesize_if(
                          node.rvalue,
                          lambda n: isinstance(n, c_ast.Assignment))
      if len(node.op)==2 and node.op.endswith("="):
        return '%s = %s %s (%s)' % (self.visit(node.lvalue), self.visit(node.lvalue), node.op[0], rval_str)
      else:
        lval_str = self.visit(node.lvalue)
        if "(void *)" in to_c(node.rvalue):
          lval_str = lval_str.replace("(","{")
          lval_str = lval_str.replace(")","}")
        return '%s %s %s' % (lval_str, node.op, rval_str)

  def visit_FuncCall(self, n):
    m = {"fabs": "abs","pmf_get_inf":"inf","pmf_sqrt":"sqrt"}
    if to_c(n.name) in m.keys():
      fref = self._parenthesize_unless_simple(n.name)
      return m[to_c(n.name)]+'(' + self.visit(n.args) + ')'
    else:
      return CommonExpressionGenerator.visit_FuncCall(self, n)
        
class SimScapeExporter(MatlabExpressionGenerator):
  def __init__(self,*args,**kwargs):
    self.cond_count = 0
    self.lvalues = []
    self.type = "double"
    self.prefix = ""
    if "cond_count" in kwargs:
      self.cond_count = kwargs["cond_count"]
      del kwargs["cond_count"]
    if "prefix" in kwargs:
      self.prefix = kwargs["prefix"]
      del kwargs["prefix"]
    if "type" in kwargs:
      self.type = kwargs["type"]
      del kwargs["type"]
    self.parent_compound = None
    MatlabExpressionGenerator.__init__(self,*args,**kwargs)

  def visit_FuncCall(self, n):
    if to_c(n.name).startswith("tlu2"):
        fref = self._parenthesize_unless_simple(n.name)
        return self.visit(n.args.exprs[0]) + " = " + "self." + fref + '(' + self.visit(n.args.exprs[1]) + ')'
    else:
      return MatlabExpressionGenerator.visit_FuncCall(self, n)

  def visit_Assignment(self, node):
    if node.lvalue.name=="out": return "% " + MatlabExpressionGenerator.visit_Assignment(self, node)
    self.lvalues.append(self.visit(node.lvalue))
    return self.prefix + MatlabExpressionGenerator.visit_Assignment(self, node)

  def visit_For(self, node):
  
       s = self.visit(node.init) + ";" + "\n"
       s += self._make_indent() + "while " + self.visit(node.cond)
       s += self._generate_stmt(node.stmt, add_indent=True)
       self.indent_level += 2
       
       assert isinstance(node.next,c_ast.UnaryOp)
       assert node.next.op=='p++'
       
       #node.stmt.show()
       s += self._make_indent() + self.visit(node.next.expr) + "=" + self.visit(node.next.expr) + "+1;\n"
       self.indent_level -= 2
       s += self._make_indent() + "end" + "\n"
       return s
       
  def visit_TernaryOp(self, n):
        s  = 'self.if_else(' + self._visit_expr(n.cond) + ','
        s += self._visit_expr(n.iftrue) + ','
        s += self._visit_expr(n.iffalse) + ')'
        return s

  def visit_If(self, node):
       #node.show()
       #ipdb.set_trace()
       
       
       #print(self.cond_count, type(self.cond_count))
       cond = "cond%d" % self.cond_count
       self.cond_count +=1
       
       # Simplify away != 0
       if isinstance(node.cond,c_ast.BinaryOp) and node.cond.op=="!=" and self.visit(node.cond.right)=="0":
          s = cond + " = " + self.visit(node.cond.left) + ";\n"
       else:
          s = cond + " = " + self.visit(node.cond) + ";\n"
       
       sp_true = SimScapeExporter(prefix=cond+"_true_" + self.prefix,cond_count=self.cond_count)
       sp_true.indent_level=self.indent_level
       st_true= sp_true.visit(node.iftrue)
       self.cond_count = sp_true.cond_count

       if node.iffalse:
         sp_false = SimScapeExporter(prefix=cond+"_false_" + self.prefix,cond_count=self.cond_count)
         sp_false.indent_level=self.indent_level
         st_false = sp_false.visit(node.iffalse)
         self.cond_count = sp_false.cond_count

       self.indent_level+=2
       if node.iffalse:
         for e in sorted(set(sp_true.lvalues) | set(sp_false.lvalues)):

           s_true = cond+"_true_"+self.prefix+e
           s_false = cond+"_false_"+self.prefix+e         
           if e not in sp_true.lvalues:
              st_true = self._make_indent() + s_true + " = " + self.prefix + e + "; % missing\n" + st_true
           if e not in sp_false.lvalues:
              st_false = self._make_indent() + s_false + " = " + self.prefix + e + "; % missing\n" + st_false
       self.indent_level-=2
       
       s+= self._make_indent() + "%"+" start block true: %s\n" % cond
       s+= st_true
       s+= self._make_indent() + "%"+" end block true: %s\n" % cond
       if node.iffalse:
         s+= self._make_indent() + "%"+" start block false: %s\n" % cond
         s+= st_false
         s+= self._make_indent() + "%"+" end block false: %s\n" % cond
         for e in sorted(set(sp_true.lvalues) | set(sp_false.lvalues)):

           s_true = cond+"_true_"+self.prefix+e
           s_false = cond+"_false_"+self.prefix+e         

           s += self._make_indent() + self.prefix+e + " = self.if_else(%s,%s,%s);\n" % (cond, s_true, s_false)
           
           self.lvalues.append(e)
       else:
         for e in sorted(set(sp_true.lvalues)):
           s_true = cond+"_true_"+ self.prefix+e
           s += self._make_indent() + self.prefix+e + " = self.if_else(%s,%s,%s);\n" % (cond, s_true, self.prefix+e)


       return s

  def visit_Return(self, n):
       return ""

  def visit_Decl(self, n, no_type=False):
       if n.name=="out": return "%" + MatlabExpressionGenerator.visit_Decl(self,n,no_type=no_type)
       if n.init is not None:
         return self.to_matlab_name(n.name) + " = [" + generator.visit(n.init) + "]"


       if self.type=="SX":
         try:
           size = n.type.dim.value
           dtype = "casadi.SX.zeros"
           try:
             if to_c(n.type.type)=="void":
               dtype = "cell"
           except:
             pass
           return self.to_matlab_name(n.name) + " = %s(%s,1)" % (dtype,size)
         except:
           pass

       return "%" + MatlabExpressionGenerator.visit_Decl(self,n,no_type=no_type)
       
  def visit_Cast(self, n):
      if self.parent_compound is not None and n in self.parent_compound.block_items:
        return "%" + MatlabExpressionGenerator.visit_Cast(self,n)
      else:
        if to_c(n).startswith("(_NeDynamicSystem"):
          return "self.extra"
        return self._parenthesize_unless_simple(n.expr)
        


  def visit_Compound(self,n):
    prev = self.parent_compound
    self.parent_compound = n
    s =  MatlabExpressionGenerator.visit_Compound(self,n)
    self.parent_compound = prev
    return s

def rvalue_to_c(node):
  return generator.visit(node)


metadata = {"variable_names": [],"input_names":[],"output_names":[],"parameter_names":[],"mmode_names":[],"n_modes":0,"s_constant_table":{},"mTable":{},"nw":0}
class MetaDataVisitor(c_ast.NodeVisitor):

  def visit_Decl(self, n, no_type=False):
       if n.name=="s_variable_data":
          for e in n.init.exprs:
            name_split = e.exprs[0].value[1:-1].split(".")
            if name_split[0] == model_file_name:
              name_split = name_split[1:] 
            #print("var", name_split, e.exprs[2].value)
            metadata["variable_names"].append(".".join(name_split))
       if n.name=="s_real_parameter_data":
          if isinstance(n.init,c_ast.ID): # Skip when NULL
            return
          for e in n.init.exprs:
            name = e.exprs[0].value[1:-1]
            metadata["parameter_names"].append(name)
       if n.name=="s_major_mode_data":
          if isinstance(n.init,c_ast.ID): # Skip when NULL
            return
          for e in n.init.exprs:
            name = e.exprs[0].value[1:-1]
            metadata["mmode_names"].append(name)
       if n.name=="s_equation_data":
          for e in n.init.exprs:
            pass #print("eq", e.exprs[0].value[1:-1].split(".")[1:], e.exprs[2].value)
            
       if n.name is not None and "s_constant_table" in n.name:
          data = []
          for e in n.init.exprs:
            data.append(float(to_c(e)))
          metadata["s_constant_table"][n.name] = data

  def visit_Assignment(self, node):
    try:
      if isinstance(node.lvalue,c_ast.StructRef):
         if node.lvalue.field.name=="mName":
           var_name = node.rvalue.value[1:-1]
           target = node.lvalue.name.name.name
           i = int(node.lvalue.name.subscript.value)
           if target=="input_info":
             assert i==len(metadata["input_names"])
             metadata["input_names"].append(var_name)
           if target=="output_info":
             assert i==len(metadata["output_names"])
             metadata["output_names"].append(var_name)
         if node.lvalue.field.name=="mNumModes":
           metadata["n_modes"]=int(node.rvalue.value)
         if node.lvalue.field.name=="mNumRealDerivedParameters":
           metadata["ndp"]=int(node.rvalue.value)
         if node.lvalue.field.name=="mNumDelays":
           metadata["nw"]=int(node.rvalue.value)
         if "mTable" in node.lvalue.field.name:
           metadata["mTable"][node.lvalue.field.name] = metadata["s_constant_table"][node.rvalue.name]
    except:
      pass
class FuncDefVisitor(c_ast.NodeVisitor):
    def visit_FuncDef(self, node):
        fields = defaultdict(list)
        name = node.decl.name
        if name.startswith(sim_file_name):
          name = name[len(sim_file_name)+1:]
        names = name.split("_")
        
        if len(names)<=1: return
        shortname = names[-2]
                
        if name.endswith("_p"):
          sp = matrices[shortname]
          if "output" in name:
            for e in node.body.block_items:
              try:
                fields[e.lvalue.field.name] = int(e.rvalue.value)
              except:
                pass
            sp.rows = fields["mNumRow"]
            sp.cols = fields["mNumCol"]
          else:
            for e in node.body.block_items:
              if not isinstance(e,c_ast.Assignment): continue
              if not isinstance(e.lvalue,c_ast.ArrayRef): continue
              fn = e.lvalue.name.field.name
              index = from_constant(e.lvalue.subscript)
              fields[fn].append((int(index), int(e.rvalue.value)))
            sp.mJc = fields["mJc"]
            sp.mIr = fields["mIr"]
            

        names_remainder = "_".join(names[1:])
        if len(names)>=2 and names[0]=="ds" and names_remainder in code_names: 
          n = names_remainder
          type = "double"
          if n in ["mode","f"]:
            type = "SX"
          generator = SimScapeExporter(type=type)
          c = codes[n] = generator.visit(node.body)
          codes_nodes[n] = node
          #node.show()
          
            
        #print('%s at %s' % (name, node.decl.coord))


code_dir = model_file_name + "_grt_rtw"

import glob

ds_file = glob.glob(os.path.join(code_dir,model_file_name+"_*_ds.c"))[0].split(os.path.sep)[1]
sim_file_name = ds_file[:-5]

ds_file_pre = os.path.join(code_dir,ds_file[:-2]+"_preprocessed.c")

with open(ds_file_pre,"w") as f_out:
  with open(os.path.join(code_dir, ds_file),"r") as f_in:
    for l in f_in.readlines():
      if l.startswith("#include") and sim_file_name+"_" in l and "ds.h" not in l:
        mod_l = l.replace(".h",".c")
        if os.path.exists(os.path.join(code_dir,mod_l.split('"')[1])):
          l = mod_l
      f_out.write(l)

cpp_path="cpp"
if os.name=='nt':
   cpp_args=["/E","/I" + basepath + "/include"]
   cpp_path=r"C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\amd64\cl.exe"
   if not os.path.exists(cpp_path):
      cpp_path=r"C:\Program Files (x86)\Microsoft Visual Studio 10.0\VC\bin\amd64\cl.exe"
else:
   cpp_args=["-I" + basepath + "/include"]

ast = parse_file(ds_file_pre, use_cpp=True, cpp_path=cpp_path,cpp_args=cpp_args)

ast = pre_transform(ast)
v = FuncDefVisitor()
v.visit(ast)


md = MetaDataVisitor()
md.visit(ast)

model_name = "Model"

constructor = []
for k in matrices:
  sp = matrices[k]

  colind = sp.colind()
  row = sp.row()
  if colind is None:
    constructor.append("self.sp_" + k + " = casadi.Sparsity({nrow},{ncol});".format(
      nrow=sp.rows,ncol=sp.cols))
  else:
    constructor.append("self.sp_" + k + " = casadi.Sparsity({nrow},{ncol},{colind},{row});".format(
      nrow=sp.rows,ncol=sp.cols,colind=colind,row=row))


standard = ["arg_x","arg_u","args_p","args_t"]

map_args = {"mode":standard+["args_q","args_w","args_s"],"f":standard+["args_q","args_w","args_s"],"update_i":["args_x","args_q","args_i"],"del_v":standard, "y": standard}
map_args_default = ["args_p"]
map_label_sys = {"arg_x": "mX","arg_u": "mU","args_p":"mP_R","args_t":"mT","args_q":"mQ","args_w":"mW","args_s":"mS"}

map_extra_body = {}

constructor.append("self.nm = "+str(metadata["n_modes"])+";")
constructor.append("self.parameter_names = {" + str(metadata["parameter_names"])[1:-1] + "};")
constructor.append("self.mmode_names = {" + str(metadata["mmode_names"])[1:-1] + "};")
constructor.append("self.variable_names = {" + str(metadata["variable_names"])[1:-1] + "};")
constructor.append("self.input_names = {" + str(metadata["input_names"])[1:-1] + "};")
constructor.append("self.output_names = {" + str(metadata["output_names"])[1:-1] + "};")
constructor.append("self.ndp = " + str(metadata["ndp"]) + ";")
constructor.append("self.nw = " + str(metadata["nw"]) + ";")
constructor.append("self.ns = 1;")

constructor.append("self.extra = struct(" + ",".join(["'"+ k+"'"+","+str(v) for k,v in metadata["mTable"].items()]) + ");")
def get_system_input_var_name(node):
  for e in node.decl.type.args.params:
    if "NeDynamicSystemInput" in e.type.type.type.names:
      return e.name

  raise Exception()

methods = []
for k in codes:
  c = codes[k]
  node = codes_nodes[k]
  
  args = map_args.get(k,map_args_default)
  extra_body = map_extra_body.get(k,[])
    
  methods.append("function ret = {name}({args})\n".format(name=k,args=",".join(["self"]+args)))
  for a in args:
    methods.append(get_system_input_var_name(node)+".%s.mX = casadi.SX(" % map_label_sys[a]+ a +");\n")
  if k in ["f"]:
    methods.append("  out.mX = casadi.SX.zeros(size(self.sp_a,1));")
  if k in ["y"]:
    methods.append("  out.mX = casadi.SX.zeros(self.ny,1);")
  if k in ["mode","f"]:
    methods.append("  out.mU = casadi.SX.zeros(size(self.sp_b,2));")
  if k in ["mode"]:
    methods.append("  out.mX = casadi.SX.zeros(self.nm,1);")
  if k in ["f"]:
    methods.append("  " + get_system_input_var_name(node)+".mM.mX = self.mode("+ ",".join(map_args["mode"])+");\n")
  if k=="del_v":
    methods.append("  out.mX = casadi.SX.zeros(self.nw,1);")
  if k=="dp_r":
    methods.append("  out.mX = casadi.SX.zeros(self.ndp,1);")
  if k!="dp_r":
    methods.append("  " + get_system_input_var_name(node) + ".mDP_R.mX = self.dp_r(args_p);")
  if k in matrices:
    methods.append("  out.mX = casadi.SX.zeros(nnz(self.sp_"+ k+ "),1);")
 
  methods+= extra_body
  methods+=c.split("\n")
  if k in matrices:
    methods.append("  ret = casadi.SX(self.sp_{name},out.mX);".format(name=k))
  else:
    methods.append("  ret = out.mX;".format(name=k))
  methods.append("end\n")

with open(model_name+".m","w") as f:
  f.write("""classdef {model_name} < DAEModel & SimscapeCasadi
    
    properties
      {properties}
      variable_names
      input_names
      output_names
      parameter_names
      mmode_names
      nm
      ndp
      extra
      nw
      ns
    end
    
    methods
      function self = {model_name}()
        self = self@DAEModel();
        {constructor}
      end
      {methods}
    end
    
end""".format(model_name=model_name,
              properties="\n      ".join("sp_" + k for k in matrices.keys()),
              constructor="\n        ".join(constructor),
              methods="\n      ".join(methods)))

# check if ode model!

#print(metadata)

