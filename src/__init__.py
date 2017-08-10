# vim: set expandtab shiftwidth=2 softtabstop=2:

# TODO - BundleAPI eh? - seems to be the class we have to override and mess with to get things imported to ChimeraX
from chimerax.core.toolshed import BundleAPI

class _MyAPI(BundleAPI):
  @staticmethod
  def get_class(class_name):
    if class_name == 'ToolUI':
      from . import tool
      return tool.ToolUI
    return None

  @staticmethod
  def start_tool(session, tool_name, **kw):
    from .tool import ToolUI
    from chimerax.core import tools
    return tools.get_singleton(session, ToolUI, 'Tempy', create=True)
 
  @staticmethod
  def register_command(command_name, logger):
    # Not sure where logger is needed or what it does yet
    from . import cmd
    # Not sure the difference between alias and register :/
    if command_name == "sccc":
      from . import cmd
      cmd.register_sccc()
    
    if command_name == "smoc":
      from . import cmd
      cmd.register_smoc()
 
    if command_name == "nmi":
      from . import cmd
      cmd.register_nmi()
    
bundle_api = _MyAPI()
