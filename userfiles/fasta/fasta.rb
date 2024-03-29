#
# CBRAIN Project
#
# Copyright (C) 2008-2021
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

# Generic model for FASTA files.
class Fasta < TextFile

  Revision_info=CbrainFileRevision[__FILE__] #:nodoc: 

  def self.pretty_type #:nodoc:
    "FASTA"
  end

  def self.file_name_pattern #:nodoc:
    /(\.fa|\.fasta)$/i
  end

  def is_viewable? #:nodoc:
    return false unless self.size.presence
    return false unless self.size < 500_000
    return false unless is_locally_synced?
    true
  end

end
