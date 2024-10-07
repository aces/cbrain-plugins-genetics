
#
# CBRAIN Project
#
# Copyright (C) 2008-2022
# The Royal Institution for the Advancement of Learning
# McGill University
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

module BoutiquesScrnaboxRessourceManager

  # Note: to access the revision info of the module,
  # you need to access the constant directly, the
  # object method revision_info() won't work.
  Revision_info=CbrainFileRevision[__FILE__] #:nodoc:

  def update_mem(step_mem, mem)
    # Update the mem according to the step_mem
    mem = step_mem > mem ? step_mem : mem
  end

  def descriptor_for_cluster_commands
    descriptor                = super.dup

    suggested_resources_by_steps = descriptor.custom_module_info("BoutiquesScrnaboxRessourceManager")

    # Steps is a string could be an integer "1" or a string "1-8"
    steps = []
    invoke_steps = invoke_params["steps"]
    if invoke_steps.include?("-")
        steps_range = invoke_steps.split("-")
        steps       = (steps_range[0].to_i..steps_range[1].to_i).to_a.map(&:to_s)
    else
        steps = [invoke_steps]
    end

    # Define walltime and mem according to data available in suggested_resources_by_steps
    walltime = 0
    mem      = 0
    for step in steps
        if suggested_resources_by_steps.has_key?(step)
            mem       = update_mem(suggested_resources_by_steps[step]["mem"], mem)
            walltime += suggested_resources_by_steps[step]["walltime"]
        end
    end

    # Add extra resources for optional flag
    if invoke_params["markergsea"]          == true && steps.start_with("7")
        mem       = update_mem(suggested_resources_by_steps[step]["mem"], mem)
        walltime += suggested_resources_by_steps["7_markergsea"]["walltime"]
    end

    if invoke_params["knownmarkers"]        == true && steps.start_with("7")
        mem       = update_mem(suggested_resources_by_steps[step]["mem"], mem)
        walltime += suggested_resources_by_steps["7_knownmarkers"]["walltime"]
    end

    if invoke_params["referenceannotation"] == true && steps.start_with("7")
        mem       = update_mem(suggested_resources_by_steps[step]["mem"], mem)
        walltime += suggested_resources_by_steps["7_referenceannotation"]["walltime"]
    end

    if invoke_params["annotate"]            == true && steps.start_with("7")
        mem       = update_mem(suggested_resources_by_steps[step]["mem"], mem)
        walltime += suggested_resources_by_steps["7_annotate"]["walltime"]
    end

    if invoke_params["addmeta"]             == true && steps.start_with("8")
        mem       = update_mem(suggested_resources_by_steps[step]["mem"], mem)
        walltime += suggested_resources_by_steps["8_addmeta"]["walltime"]
    end

    if invoke_params["rundge"]              == true && steps.start_with("8")
        mem       = update_mem(suggested_resources_by_steps[step]["mem"], mem)
        walltime += suggested_resources_by_steps["8_rundge"]["walltime"]
    end

    # In minute in descriptor, need to be converted in second
    walltime = walltime * 60
    descriptor["suggested-resources"]["walltime-estimate"] = walltime
    descriptor["suggested-resources"]["mem"]               = mem

    return descriptor
  end

end