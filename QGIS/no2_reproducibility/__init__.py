def classFactory(iface):
    from .plugin import NO2Plugin
    return NO2Plugin(iface)
