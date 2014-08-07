# vim: set expandtab ts=4 sw=4:

class CmdLine:
    
    SIZE = (500, 25)

    def __init__(self):
        import wx
        from .tool_api import ToolWindow
        self.tool_window = ToolWindow("Command Line", "General",
            size=self.SIZE, placement=wx.BOTTOM)
        parent = self.tool_window.ui_area
        self.text = wx.TextCtrl(parent, size=self.SIZE,
            style=wx.TE_PROCESS_ENTER | wx.TE_NOHIDESEL)
        self.text.Bind(wx.EVT_TEXT_ENTER, self.OnEnter)

    def OnEnter(self, event):
        cmd = self.text.GetLineText(0)
        self.text.SelectAll()
