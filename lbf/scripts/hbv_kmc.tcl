#!/usr/bin/tclsh

#Render hbv with variable atoms and bonds
#requires topotools

set env(VMDNODISRUPTHACK) 1 

proc load_traj {myfile} {
  #Loop through trajlammps.dat and extract
  #each configuration as a separate molecule

  #Try opening file
  if { [catch {open $myfile r} fid] } {
    puts stderr "Trajectory file does not exist.\n"
    #exit 1
  }
  set contents [read -nonewline $fid]
  close $fid
  set splitCont [split $contents "\n"] 

  set linecount 0
  
  #Get rid of temp file if it exists
  if { [file exists "temp.dat"] == 1} {
    file delete [open "temp.dat"]
  }

  #Create temp file and write out each configuration
  #in the trajectory to it
  set tempfile [open "temp.dat" a]
  foreach ele $splitCont {
    #The "LAMMPS" substring denotes the header of
    #a new configuration
    if {[string first "LAMMPS" $ele] != -1} {
      if { $linecount != 0} {
        close $tempfile
        #Load the configuration into a new molecule
        #topo readlammpsdata "temp.dat"
        load_conf "temp.dat"
      }
      set linecount 0
      puts "$ele"
      close $tempfile
      file delete temp.dat
      set tempfile [open "temp.dat" w]
    }
    puts $tempfile "$ele\n"
    incr linecount
  }
  file delete temp.dat

  set index 0
  foreach id [molinfo list] {
    if { $index != 0} {
      mol off $id
    }
    incr index
  }
  set index 0
}

proc load_conf {myfile} {

  if { [catch {open $myfile r} fid] } {
    puts stderr "File does not exist.\n"
    #exit 1
  }
  topo readlammpsdata $myfile

  mol delrep 0 top

  #A conformation (red)
  mol color ColorID 0
  mol representation Licorice 0.12 12 12
  mol selection name 2
  mol addrep top

  #B conformation (white)
  mol color ColorID 8
  mol representation Licorice 0.12 12 12
  mol selection name 3
  mol addrep top

  #C conformation (cyan)
  mol color ColorID 10
  mol representation Licorice 0.12 12 12
  mol selection name 1
  mol addrep top

  #D conformation (dark turquoise)
  mol color ColorID 0
  mol representation Licorice 0.12 12 12
  mol selection name 4
  mol addrep top

  #Drug molecule
  mol color ColorID 3
  mol representation VDW 0.1 12.000000
  mol material Opaque
  mol selection name 5
  mol addrep top

  color Display Background white
  axes location off
  display projection Orthographic
  #scale to 0.04

}

proc make_movie {path interv} {

  set index 0
  foreach id [molinfo list] {
    mol off $id
  }
  display update on
  foreach id [molinfo list] {
    if { [expr $id - [lindex [molinfo list] 0]] % $interv == 0 } {
      mol on $id
      mol top $id
      display resetview
      scale to 0.4
      display update

      file mkdir $path ;# [file mkdir] in Tcl is like mkdir -p
      file mkdir "$path/images"
      render TachyonInternal [format "$path/images/image_%d.tga" $index]
      incr index
      mol off $id
    }
  }

}

proc myanimate {waittime} {
  foreach id [molinfo list] {
    mol off $id
  }
  display update on
  #wait for 1 second
  after 1000 
  set index 0
  foreach id [molinfo list] {
    puts "$index"
    #if { $index != 0 } {
    #  puts "[$id - 1]"
    #  mol off [lindex [molinfo list] $index]
    #}
    mol on $id
    mol top $id
    display resetview
    scale to 0.4
    display update
    after $waittime
    mol off $id
    incr index
  }
}

