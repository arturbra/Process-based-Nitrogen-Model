Pesz_no3 = fPesz(UFi_sz, D_no3)

if m_usz < (n - 1):
    if Pesz_no3 <= 2:
        if sz_layer == 0:  # first cell
            dc_no3 = cjplus1_no3 - 2 * concentration_sz + concentration_last_layer_usz
            dc_dz_no3 = (cjplus1_no3 - concentration_last_layer_usz) / (2 * dz)

        elif sz_layer == (m_sz - 1):  # last cell
            dc_no3 = concentration_sz - 2 * concentration_sz + concentration_sz_layer_above
            dc_dz_no3 = (concentration_sz - concentration_sz_layer_above) / dz

        else:
            dc_no3 = cjplus1_no3 - 2 * concentration_sz + concentration_sz_layer_above
            dc_dz_no3 = (cjplus1_no3 - concentration_sz_layer_above) / (2 * dz)

    else:  # Pesz > 2
        if sz_layer == 0:  # first cell
            dc_no3 = cjplus1_no3 - 2 * concentration_sz + concentration_last_layer_usz
            dc_dz_no3 = (concentration_sz - concentration_last_layer_usz) / dz

        elif sz_layer == (m_sz - 1):  # last cell
            dc_no3 = concentration_sz - 2 * concentration_sz + concentration_sz_layer_above
            dc_dz_no3 = (concentration_sz - concentration_sz_layer_above) / dz

        else:
            dc_no3 = cjplus1_no3 - 2 * concentration_sz + concentration_sz_layer_above
            dc_dz_no3 = (concentration_sz - concentration_sz_layer_above) / dz

if m_usz == (n - 1):
    dc_no3 = concentration_sz - 2 * concentration_sz + concentration_sz_layer_above
    dc_dz_no3 = (concentration_sz - concentration_sz_layer_above) / (2 * dz)

if m_usz == 0:
    if Pesz_no3 <= 2:
        if sz_layer == 0:  # first cell
            dc_no3 = cjplus1_no3 - 2 * concentration_sz + cpi_no3
            dc_dz_no3 = (cjplus1_no3 - cpi_no3) / (2 * dz)

        elif sz_layer == (m_sz - 1):  # last cell
            dc_no3 = concentration_sz - 2 * concentration_sz + concentration_sz_layer_above
            dc_dz_no3 = (concentration_sz - concentration_sz_layer_above) / dz

        else:
            dc_no3 = cjplus1_no3 - 2 * concentration_sz + concentration_sz_layer_above
            dc_dz_no3 = (cjplus1_no3 - concentration_sz_layer_above) / (2 * dz)

    else:  # Pusz > 2
        if sz_layer == 0:  # first cell
            dc_no3 = cjplus1_no3 - 2 * concentration_sz + cpi_no3
            dc_dz_no3 = (concentration_sz - cpi_no3) / dz

        elif sz_layer == (m_sz - 1):  # last cell
            dc_no3 = concentration_sz - 2 * concentration_sz + concentration_sz_layer_above
            dc_dz_no3 = (concentration_sz - concentration_sz_layer_above) / dz

        else:
            dc_no3 = cjplus1_no3 - 2 * concentration_sz + concentration_sz_layer_above
            dc_dz_no3 = (concentration_sz - concentration_sz_layer_above) / dz

delta_c_no3 = ftransp(teta_b_i, teta_b_iplus1, concentration_sz, cs_no3, dc_no3, dc_dz_no3, 0, 0, D_no3, UFi_sz, Rxi_3_no3)
if teta_b_iplus1 > 0:
    ci1_no3 = concentration_sz + delta_c_no3
else:
    ci1_no3 = 0

if ci1_no3 <= 0.0000000000000001:
    ci1_no3 = 0
else:
    ci1_no3 = ci1_no3

cj_i1_no3.append(ci1_no3)
