from qgis.core import QgsApplication
from .provider import NO2Provider
class NO2Plugin:
    def __init__(self,iface):self.iface=iface;self.provider=None
    def initGui(self):
        self.provider=NO2Provider();QgsApplication.processingRegistry().addProvider(self.provider)
    def unload(self):
        if self.provider:QgsApplication.processingRegistry().removeProvider(self.provider);self.provider=None
