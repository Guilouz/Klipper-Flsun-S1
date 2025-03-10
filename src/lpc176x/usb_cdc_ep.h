#ifndef __LPC176X_USB_CDC_EP_H
#define __LPC176X_USB_CDC_EP_H

enum {
    USB_CDC_EP_ACM = 1,
    USB_CDC_EP_BULK_OUT = 2,
    USB_CDC_EP_BULK_IN = 5,
};

enum {
    USB_CDC_EP0_SIZE = 16,
    USB_CDC_EP_ACM_SIZE = 8,
    USB_CDC_EP_BULK_OUT_SIZE = 64,
    USB_CDC_EP_BULK_IN_SIZE = 64,
};

#endif // usb_cdc_ep.h
