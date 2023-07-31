import math


def get_pred(Wp, Wm, W, Z, corr):

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

    WZ = {}
    WZ['val'] = W['val']/Z['val']
    WZ['err_scale_up'] = WZ['val'] * math.sqrt( (W['err_scale_up']/W['val'])**2 + (Z['err_scale_up']/Z['val'])**2 - 2*corr['WZ_scale']*W['err_scale_up']*Z['err_scale_up']/(W['val']*Z['val']))
    WZ['err_scale_down'] = WZ['val'] * math.sqrt( (W['err_scale_down']/W['val'])**2 + (Z['err_scale_down']/Z['val'])**2 - 2*corr['WZ_scale']*W['err_scale_down']*Z['err_scale_down']/(W['val']*Z['val']))
    WZ['err_pdf_up'] = WZ['val'] * math.sqrt( (W['err_pdf_up']/W['val'])**2 + (Z['err_pdf_up']/Z['val'])**2 - 2*corr['WZ_pdf']*W['err_pdf_up']*Z['err_pdf_up']/(W['val']*Z['val']))
    WZ['err_pdf_down'] = WZ['val'] * math.sqrt( (W['err_pdf_down']/W['val'])**2 + (Z['err_pdf_down']/Z['val'])**2 - 2*corr['WZ_pdf']*W['err_pdf_down']*Z['err_pdf_down']/(W['val']*Z['val']))

    print()

    # print(Wp)
    # print(Wm)
    # print(W)
    # print(Z)

    # print()

    # print(WpWm)
    # print(WpZ)
    # print(WmZ)
    # print(WZ)

    # print()

    print(f"Wp: {Wp['val']:2.1f} +{Wp['err_pdf_up']:2.1f} -{Wp['err_pdf_down']:2.1f} (PDF) +{Wp['err_scale_up']:2.1f} -{Wp['err_scale_down']:2.1f} (scale)")
    print(f"Wm: {Wm['val']:2.1f} +{Wm['err_pdf_up']:2.1f} -{Wm['err_pdf_down']:2.1f} (PDF) +{Wm['err_scale_up']:2.1f} -{Wm['err_scale_down']:2.1f} (scale)")
    print(f"W: {W['val']:2.1f} +{W['err_pdf_up']:2.1f} -{W['err_pdf_down']:2.1f} (PDF) +{W['err_scale_up']:2.1f} -{W['err_scale_down']:2.1f} (scale)")
    print(f"Z: {Z['val']:2.1f} +{Z['err_pdf_up']:2.1f} -{Z['err_pdf_down']:2.1f} (PDF) +{Z['err_scale_up']:2.1f} -{Z['err_scale_down']:2.1f} (scale)")

    print()

    print(f"WpZ: {WpZ['val']:2.3f} +{WpZ['err_pdf_up']:2.3f} -{WpZ['err_pdf_down']:2.3f} (PDF) +{WpZ['err_scale_up']:2.3f} -{WpZ['err_scale_down']:2.3f} (scale)")
    print(f"WmZ: {WmZ['val']:2.3f} +{WmZ['err_pdf_up']:2.3f} -{WmZ['err_pdf_down']:2.3f} (PDF) +{WmZ['err_scale_up']:2.3f} -{WmZ['err_scale_down']:2.3f} (scale)")
    print(f"WpWm: {WpWm['val']:2.3f} +{WpWm['err_pdf_up']:2.3f} -{WpWm['err_pdf_down']:2.3f} (PDF) +{WpWm['err_scale_up']:2.3f} -{WpWm['err_scale_down']:2.3f} (scale)")
    print(f"WZ: {WZ['val']:2.3f} +{WZ['err_pdf_up']:2.3f} -{WZ['err_pdf_down']:2.3f} (PDF) +{WZ['err_scale_up']:2.3f} -{WZ['err_scale_down']:2.3f} (scale)")

    print()


### fid
Wp_fid = {}
Wm_fid = {}
W_fid = {}
Z_fid = {}

Wp_fid['val'] = 5427.8
Wm_fid['val'] = 4155.3
W_fid['val'] = Wp_fid['val'] + Wm_fid['val']
Z_fid['val'] = 753.05

Wp_fid['err_scale_up'] = 0
Wm_fid['err_scale_up'] = 0
W_fid['err_scale_up'] = Wp_fid['err_scale_up'] + Wm_fid['err_scale_up']
Z_fid['err_scale_up'] = 2.07

Wp_fid['err_scale_down'] = 0
Wm_fid['err_scale_down'] = 0
W_fid['err_scale_down'] = Wp_fid['err_scale_down'] + Wm_fid['err_scale_down']
Z_fid['err_scale_down'] = 4.47

Wp_fid['err_pdf_up'] = 77.80
Wm_fid['err_pdf_up'] = 68.3
W_fid['err_pdf_up'] = Wp_fid['err_pdf_up'] + Wm_fid['err_pdf_up']
Z_fid['err_pdf_up'] = 14.1

Wp_fid['err_pdf_down'] = 121.9
Wm_fid['err_pdf_down'] = 107.8
W_fid['err_pdf_down'] = Wp_fid['err_pdf_down'] + Wm_fid['err_pdf_down']
Z_fid['err_pdf_down'] = 22.6

corr_fid = {}
corr_fid['WW_scale'] = 1
corr_fid['WZ_scale'] = 0
corr_fid['WW_pdf'] = 1
corr_fid['WZ_pdf'] = 0


### total
Wp_tot = {}
Wm_tot = {}
W_tot = {}
Z_tot = {}

Wp_tot['val'] = 11939.2
Wm_tot['val'] = 8857.4
W_tot['val'] = Wp_tot['val'] + Wm_tot['val']
Z_tot['val'] = 2019.81

Wp_tot['err_scale_up'] = 0
Wm_tot['err_scale_up'] = 0
W_tot['err_scale_up'] = Wp_tot['err_scale_up'] + Wm_tot['err_scale_up']
Z_tot['err_scale_up'] = 17.6

Wp_tot['err_scale_down'] = 0
Wm_tot['err_scale_down'] = 0
W_tot['err_scale_down'] = Wp_tot['err_scale_down'] + Wm_tot['err_scale_down']
Z_tot['err_scale_down'] = 23.3

Wp_tot['err_pdf_up'] = 268.0
Wm_tot['err_pdf_up'] = 184.3
W_tot['err_pdf_up'] = Wp_tot['err_pdf_up'] + Wm_tot['err_pdf_up']
Z_tot['err_pdf_up'] = 34.0

Wp_tot['err_pdf_down'] = 416.2
Wm_tot['err_pdf_down'] = 227.6
W_tot['err_pdf_down'] = Wp_tot['err_pdf_down'] + Wm_tot['err_pdf_down']
Z_tot['err_pdf_down'] = 50.6

corr_tot = {}
corr_tot['WW_scale'] = 1
corr_tot['WZ_scale'] = 0
corr_tot['WW_pdf'] = 1
corr_tot['WZ_pdf'] = 0


print('fiducial pred:')
get_pred(Wp_fid, Wm_fid, W_fid, Z_fid, corr_fid)

print('total pred:')
get_pred(Wp_tot, Wm_tot, W_tot, Z_tot, corr_tot)
