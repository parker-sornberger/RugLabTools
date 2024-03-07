draw delete all
set filename modesvmd.nmm
set data [open $filename r]
set vmd_pick_mol 0
draw color red
set scale 10

foreach line [split [read $data] \n] {
		if {[lindex $line 0] == "ATOM"} {
			set vmd_pick_atom [expr {[lindex $line 1] - 1}]

			set sel [atomselect $vmd_pick_mol "index $vmd_pick_atom"]
			set coords1 [lindex [$sel get {x y z}] 0]

			set coords3 "[lindex $line 2] [lindex $line 3] [lindex $line 4]"
			set coords3 [vecscale $scale $coords3]
			set coords2 [vecscale 0.6 $coords3]
			set coords2 [vecadd $coords1 $coords2]
			set coords3 [vecadd $coords1 $coords3]
	
			draw cylinder $coords1 $coords2 radius 0.075
			draw cone $coords2 $coords3 radius 0.2

			puts stdout "Atom $vmd_pick_atom on Molecule $vmd_pick_mol done."
			puts stdout "  VECTOR: $coords1 to $coords2 to $coords3"
		} 
	}	
	close $data
	
