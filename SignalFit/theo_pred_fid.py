import math

Wp = {}
Wm = {}
Z = {}

Wp['val'] = 5427.8
Wm['val'] = 4155.3
Z['val'] = 753.05

Wp['err_scale_up'] = 0
Wm['err_scale_up'] = 0
Z['err_scale_up'] = 2.07

Wp['err_scale_down'] = 0
Wm['err_scale_down'] = 0
Z['err_scale_down'] = 4.47

Wp['err_pdf_up'] = 77.80
Wm['err_pdf_up'] = 68.3
Z['err_pdf_up'] = 14.1

Wp['err_pdf_down'] = 121.9
Wm['err_pdf_down'] = 107.8
Z['err_pdf_down'] = 22.6

corr = {}
corr['WW_scale'] = 1
corr['WZ_scale'] = 0
corr['WW_pdf'] = 1
corr['WZ_pdf'] = 0


WpWm = {}
WpWm['val'] = Wp['val']/Wm['val']
WpWm['err_scale_up'] = WpWm['val'] * math.sqrt( (Wp['err_scale_up']/Wp['val'])**2 + (Wm['err_scale_up']/Wm['val'])**2 - 2*corr['WW_scale']*Wp['err_scale_up']*Wm['err_scale_up']/(Wp['val']*Wm['val']))
WpWm['err_scale_down'] = WpWm['val'] * math.sqrt( (Wp['err_scale_down']/Wp['val'])**2 + (Wm['err_scale_down']/Wm['val'])**2 - 2*corr['WW_scale']*Wp['err_scale_down']*Wm['err_scale_down']/(Wp['val']*Wm['val']))
WpWm['err_pdf_up'] = WpWm['val'] * math.sqrt( (Wp['err_pdf_up']/Wp['val'])**2 + (Wm['err_pdf_up']/Wm['val'])**2 - 2*corr['WW_pdf']*Wp['err_pdf_up']*Wm['err_pdf_up']/(Wp['val']*Wm['val']))
WpWm['err_pdf_down'] = WpWm['val'] * math.sqrt( (Wp['err_pdf_down']/Wp['val'])**2 + (Wm['err_pdf_down']/Wm['val'])**2 - 2*corr['WW_pdf']*Wp['err_pdf_down']*Wm['err_pdf_down']/(Wp['val']*Wm['val']))

WpZ = {}
WpZ['val'] = Wp['val']/Z['val']
WpZ['err_scale_up'] = WpZ['val'] * math.sqrt( (Wp['err_scale_up']/Wp['val'])**2 + (Z['err_scale_up']/Z['val'])**2 - 2*corr['WZ_scale']*Wp['err_scale_up']*Z['err_scale_up']/(Wp['val']*Z['val']))
WpZ['err_scale_down'] = WpZ['val'] * math.sqrt( (Wp['err_scale_down']/Wp['val'])**2 + (Z['err_scale_down']/Z['val'])**2 - 2*corr['WZ_scale']*Wp['err_scale_down']*Z['err_scale_down']/(Wp['val']*Z['val']))
WpZ['err_pdf_up'] = WpZ['val'] * math.sqrt( (Wp['err_pdf_up']/Wp['val'])**2 + (Z['err_pdf_up']/Z['val'])**2 - 2*corr['WZ_pdf']*Wp['err_pdf_up']*Z['err_pdf_up']/(Wp['val']*Z['val']))
WpZ['err_pdf_down'] = WpZ['val'] * math.sqrt( (Wp['err_pdf_down']/Wp['val'])**2 + (Z['err_pdf_down']/Z['val'])**2 - 2*corr['WZ_pdf']*Wp['err_pdf_down']*Z['err_pdf_down']/(Wp['val']*Z['val']))

WmZ = {}
WmZ['val'] = Wm['val']/Z['val']
WmZ['err_scale_up'] = WmZ['val'] * math.sqrt( (Wm['err_scale_up']/Wm['val'])**2 + (Z['err_scale_up']/Z['val'])**2 - 2*corr['WZ_scale']*Wm['err_scale_up']*Z['err_scale_up']/(Wm['val']*Z['val']))
WmZ['err_scale_down'] = WmZ['val'] * math.sqrt( (Wm['err_scale_down']/Wm['val'])**2 + (Z['err_scale_down']/Z['val'])**2 - 2*corr['WZ_scale']*Wm['err_scale_down']*Z['err_scale_down']/(Wm['val']*Z['val']))
WmZ['err_pdf_up'] = WmZ['val'] * math.sqrt( (Wm['err_pdf_up']/Wm['val'])**2 + (Z['err_pdf_up']/Z['val'])**2 - 2*corr['WZ_pdf']*Wm['err_pdf_up']*Z['err_pdf_up']/(Wm['val']*Z['val']))
WmZ['err_pdf_down'] = WmZ['val'] * math.sqrt( (Wm['err_pdf_down']/Wm['val'])**2 + (Z['err_pdf_down']/Z['val'])**2 - 2*corr['WZ_pdf']*Wm['err_pdf_down']*Z['err_pdf_down']/(Wm['val']*Z['val']))

print()

print(Wp)
print(Wm)
print(Z)

print()

print(WpWm)
print(WpZ)
print(WmZ)

print()

print(f"Wp: {Wp['val']:2.1f} +{Wp['err_pdf_up']:2.1f} -{Wp['err_pdf_down']:2.1f} (PDF) +{Wp['err_scale_up']:2.1f} -{Wp['err_scale_down']:2.1f} (scale)")
print(f"Wm: {Wm['val']:2.1f} +{Wm['err_pdf_up']:2.1f} -{Wm['err_pdf_down']:2.1f} (PDF) +{Wm['err_scale_up']:2.1f} -{Wm['err_scale_down']:2.1f} (scale)")
print(f"Z: {Z['val']:2.1f} +{Z['err_pdf_up']:2.1f} -{Z['err_pdf_down']:2.1f} (PDF) +{Z['err_scale_up']:2.1f} -{Z['err_scale_down']:2.1f} (scale)")

print()

print(f"WpZ: {WpZ['val']:2.3f} +{WpZ['err_pdf_up']:2.3f} -{WpZ['err_pdf_down']:2.3f} (PDF) +{WpZ['err_scale_up']:2.3f} -{WpZ['err_scale_down']:2.3f} (scale)")
print(f"WmZ: {WmZ['val']:2.3f} +{WmZ['err_pdf_up']:2.3f} -{WmZ['err_pdf_down']:2.3f} (PDF) +{WmZ['err_scale_up']:2.3f} -{WmZ['err_scale_down']:2.3f} (scale)")
print(f"WpWm: {WpWm['val']:2.3f} +{WpWm['err_pdf_up']:2.3f} -{WpWm['err_pdf_down']:2.3f} (PDF) +{WpWm['err_scale_up']:2.3f} -{WpWm['err_scale_down']:2.3f} (scale)")
