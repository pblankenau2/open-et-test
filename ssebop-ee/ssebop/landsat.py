import ee


system_properties = ['system:index', 'system:time_start']

def prep_func(image):
    """Compute LST and NDVI for a prepped Landsat image

    This function can be mapped over a collection

    Input EE Image must have the following bands:
        'nir', 'red', 'thermal'
    Input EE Image must have the following joined properties/images:
        'k1_constant'/'k2_constant' (from Landsat prep functions)
    """
    ndvi = ndvi_func(image)
    lst = lst_func(image)
    return ee.Image([lst, ndvi]).rename(['lst', 'ndvi']) \
        .copyProperties(image, system_properties)


def ndvi_func(image):
    """Compute NDVI for a prepped image

    This function can be mapped over a collection

    Input EE Image must have the following bands:
        'nir' and 'red'
    """
    return ee.Image(image).normalizedDifference(['nir', 'red']) \
        .rename(['ndvi'])


def lst_func(image):
    """Compute LST for a prepped image

    This function can be mapped over a collection

    Input EE Image must have the following bands:
        'nir', 'red', 'thermal'
    Input EE Image must have the following joined properties/images:
        'k1_constant'/'k2_constant' (from Landsat prep functions)
    """
    # Get properties from image
    ts_brightness = ee.Image(image).select(['thermal'])
    emissivity = _emissivity(
        ee.Image(image).normalizedDifference(['nir', 'red']))
    k1 = ee.Number(image.get('k1_constant'))
    k2 = ee.Number(image.get('k2_constant'))

    # First back out radiance from brightness temperature
    # Then recalculate emissivity corrected Ts
    thermal_rad_toa = ts_brightness.expression(
        'k1 / (exp(k2 / ts_brightness) - 1)',
        {'ts_brightness': ts_brightness, 'k1': k1, 'k2': k2})

    # tnb = 0.866   # narrow band transmissivity of air
    # rp = 0.91     # path radiance
    # rsky = 1.32   # narrow band clear sky downward thermal radiation
    rc = thermal_rad_toa.expression(
        '((thermal_rad_toa - rp) / tnb) - ((1. - emiss) * rsky)',
        {
            'thermal_rad_toa': thermal_rad_toa,
            'emiss': emissivity,
            'rp': 0.91, 'tnb': 0.866, 'rsky': 1.32})
    lst = rc.expression(
        'k2 / log(emiss * k1 / rc + 1)',
        {'emiss': emissivity, 'rc': rc, 'k1': k1, 'k2': k2})

    return lst.rename(['lst'])


def _emissivity(ndvi):
    """Compute emissivity from NDVI

    This function can be mapped over a collection
    """
    Pv = ndvi.expression(
        '((ndvi - 0.2) / 0.3) ** 2', {'ndvi': ndvi})
    # ndviRangevalue = ndvi_image.where(
    #     ndvi_image.gte(0.2).And(ndvi_image.lte(0.5)), ndvi_image)
    # Pv = ndviRangevalue.expression(
    # '(((ndviRangevalue - 0.2)/0.3)**2',{'ndviRangevalue':ndviRangevalue})

    # Assuming typical Soil Emissivity of 0.97 and Veg Emissivity of 0.99
    #   and shape Factor mean value of 0.553
    dE = Pv.expression(
        '(((1 - 0.97) * (1-Pv)) * (0.55 * 0.99))', {'Pv': Pv})
    RangeEmiss = dE.expression(
        '((0.99 * Pv) + (0.97 * (1 - Pv)) + dE)', {'Pv': Pv, 'dE': dE})

    # RangeEmiss = 0.989 # dE.expression(
    #  '((0.99*Pv)+(0.97 *(1-Pv))+dE)',{'Pv':Pv, 'dE':dE})
    emissivity = ndvi \
        .where(ndvi.lt(0), 0.985) \
        .where((ndvi.gte(0)).And(ndvi.lt(0.2)), 0.977) \
        .where(ndvi.gt(0.5), 0.99) \
        .where((ndvi.gte(0.2)).And(ndvi.lte(0.5)), RangeEmiss)
    emissivity = emissivity.clamp(0.977, 0.99)
    return emissivity.select([0], ['emissivity']) \
        .copyProperties(ndvi, system_properties)
