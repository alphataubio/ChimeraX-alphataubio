# vim: set expandtab ts=4 sw=4:

class ToolWindow:

    placements = ["right", "left", "top", "bottom"]

    def __init__(self, toolName, category, session, menus=False,
            prefer_detached=False, icon=None, size=None, placement=None,
            destroy_hides=False):
        try:
            self.__toolkit = _Wx(self, toolName, menus, session,
                prefer_detached, size, placement, destroy_hides)
        except ImportError:
            # browser version
            raise NotImplementedError("Browser tool API not implemented")
        self.ui_area = self.__toolkit.ui_area

    def destroy(self):
        self.__toolkit.destroy()
        self.__toolkit = None

    def getShown(self):
        return self.__toolkit.shown

    def setShown(self, shown):
        self.__toolkit.shown = shown

    shown = property(getShown, setShown)

class _Wx:

    def __init__(self, tool_window, toolName, menus, session, prefer_detached,
            size, placement, destroy_hides):
        self.tool_window = tool_window
        self.destroy_hides = destroy_hides
        import wx
        self.main_window = mw = session.main_window
        if not mw:
            raise RuntimeError("No main window or main window dead")
        if size is None:
            size = wx.DefaultSize
        class WxToolPanel(wx.Panel):
            def __init__(self, parent, destroy_hides=destroy_hides, **kw):
                self._destroy_hides = destroy_hides
                wx.Panel.__init__(self, parent, **kw)
                self.Bind(wx.EVT_CLOSE, self.OnClose)

            def OnClose(self, event):
                if self.destroy_hides and event.CanVeto():
                    self.Shown = False
                else:
                    self.Destroy()

        self.ui_area = WxToolPanel(mw, name=toolName, size=size)
        placements = self.tool_window.placements
        if placement is None:
            placement = wx.RIGHT
        elif placement not in placements:
            raise ValueError("placement value must be one of: {}".format(
                ", ".join(placements)))
        else:
            placement = dict(zip(placements, [wx.RIGHT, wx.LEFT,
                wx.TOP, wx.BOTTOM]))[placement]
        mw.aui_mgr.AddPane(self.ui_area, placement, toolName)
        mw.aui_mgr.Update()
        if prefer_detached:
           mw.aui_mgr.GetPane(self.ui_area).Float()
        self._pane_info = None

    def destroy(self):
        if self.tool_window.ui_panel:
            self.shown = False
            if self.destroy_hides:
                return
            self.ui_area.Destroy()
        self.tool_window = None

    def getShown(self):
        return self.ui_area.Shown

    def setShown(self, shown):
        aui_mgr = self.main_window.aui_mgr
        if shown:
            if self._pane_info:
                # has been hidden at least once
                aui_mgr.AddPane(self.ui_area, self._pane_info)
                self._pane_info = None
        else:
            self._pane_info = aui_mgr.GetPane(self.ui_area)
            aui_mgr.DetachPane(self.ui_area)

        self.ui_area.Shown = shown

    shown = property(getShown, setShown)
