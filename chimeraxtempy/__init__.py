# vim: set expandtab shiftwidth=2 softtabstop=2:

from chimerax.core.toolshed import BundleAPI

class _MyAPI(BundleAPI):
  @staticmethod
  def get_class(class_name):
    #if class_name == 'ToolUI':
    #  from . import tool
    #  return tool.ToolUI
    return None

  @staticmethod
  def start_tool(session, tool_name, **kw):
    #from .tool import ToolUI
    #return ToolUI(session, tool_name, **kw)
    pass

  @staticmethod
  def register_command(command_name):
    from . import cmd
    # Not sure the difference between alias and register :/
    if command_name == "tempy":
      from chimerax.core.commands import register
      register(command_name, cmd.tempy_desc, cmd.tempy)
    elif command_name == "tempytest":
      from chimerax.core.commands import create_alias
      create_alias(command_name, "tempytest $*")

bundle_api = _MyAPI()

def main():
  """Entry point for the application script"""
  print("Call your main application code here")
