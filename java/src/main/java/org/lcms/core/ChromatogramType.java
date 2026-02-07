package org.lcms.core;

/**
 * Type of chromatogram.
 */
public enum ChromatogramType {
    UNKNOWN,
    TIC,    // Total Ion Current
    BPC,    // Base Peak Chromatogram
    XIC,    // Extracted Ion Chromatogram
    SRM,    // Selected Reaction Monitoring
    MRM,    // Multiple Reaction Monitoring
    SIM     // Selected Ion Monitoring
}
