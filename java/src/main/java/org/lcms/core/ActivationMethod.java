package org.lcms.core;

/**
 * MS/MS activation methods.
 */
public enum ActivationMethod {
    UNKNOWN,
    CID,    // Collision-Induced Dissociation
    HCD,    // Higher-energy Collisional Dissociation
    ETD,    // Electron Transfer Dissociation
    ECD,    // Electron Capture Dissociation
    UVPD,   // Ultraviolet Photodissociation
    IRMPD   // Infrared Multiphoton Dissociation
}
