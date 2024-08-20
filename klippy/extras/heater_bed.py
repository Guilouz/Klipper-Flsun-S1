# Support for a heated bed
#
# Copyright (C) 2018-2019  Kevin O'Connor <kevin@koconnor.net>
#
# This file may be distributed under the terms of the GNU GPLv3 license.

class PrinterHeaterBed:
    def __init__(self, config):
        self.printer = config.get_printer()
        pheaters = self.printer.load_object(config, 'heaters')
        self.heater = pheaters.setup_heater(config, 'B')
        self.get_status = self.heater.get_status
        self.stats = self.heater.stats
        # Register commands
        gcode = self.printer.lookup_object('gcode')
        gcode.register_command("M140", self.cmd_M140)
        gcode.register_command("M190", self.cmd_M190)
    def cmd_M140(self, gcmd, wait=False):
        # Set Bed Temperature
        temp = gcmd.get_float('S', 0.)
        # Start FLSUN Changes
        hotbed_0 = gcmd.get_float('A', 1)
        hotbed_1 = gcmd.get_float('B', 1)
        # End FLSUN Changes
        pheaters = self.printer.lookup_object('heaters')
        # Start FLSUN Changes
        #pheaters.set_temperature(self.heater, temp, wait)
        gcode = self.printer.lookup_object('gcode')
        if wait:
            if hotbed_0 == 0:
                pheaters.set_temperature(self.heater, temp, False)
            if hotbed_1 == 0:
                gcode.run_script_from_command("SET_HEATER_TEMPERATURE HEATER=heater_bed_2 TARGET=%f WAIT=0" % temp)
        if hotbed_0==1:
            pheaters.set_temperature(self.heater, temp, wait)
        if hotbed_1==1:
            if wait:
                gcode.run_script_from_command("SET_HEATER_TEMPERATURE HEATER=heater_bed_2 TARGET=%f WAIT=1" % temp) 
            else:
                gcode.run_script_from_command("SET_HEATER_TEMPERATURE HEATER=heater_bed_2 TARGET=%f WAIT=0" % temp)
        if wait:
            if hotbed_0 == 0:
                pheaters.set_temperature(self.heater, 0, False)
            if hotbed_1 == 0:
                gcode.run_script_from_command("SET_HEATER_TEMPERATURE HEATER=heater_bed_2 TARGET=%f WAIT=0" % 0)
        # End FLSUN Changes
    def cmd_M190(self, gcmd):
        # Set Bed Temperature and Wait
        self.cmd_M140(gcmd, wait=True)

def load_config(config):
    return PrinterHeaterBed(config)
