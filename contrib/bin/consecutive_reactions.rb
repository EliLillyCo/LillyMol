#!/usr/bin/env ruby

# Sometimes we need to perform the same reaction N times. For
# example with positional analogue scanning, we may wish
# to attach precisely 2 functional groups to a molecule.
# I am not sure this can be done with a single invocation of trxn.


require_relative('lib/iwcmdline')

def usage
  msg = <<-MSG
Runs consecutive trxn invocations, so that a given reaction can be applied
to the same molecule N times. Certain types of Positional Analogue Scanning.

The first trxn invocation might add a given atom to the molecule in one location.
The second invocation then adds another atom, at a different location.
If this script is run with `-N 2` there will be two pipelined trxn invocations and the
resulting molecules will have two functional groups added.

For example if you with to add two fluorines to aroamtic carbon sites in a molecule start
with a reaction which adds a single F atom like (id's omitted)
```
scaffold {
  smarts: "cH"
}
sidechain {
  reagent: "F"
  smarts: "F"
  join { a1: 0 a2: 0}
}
```
The invocation might look like
```
consecutive_reactions -P file.rxn -N 2 file.smi > new.smi
```
If instead you wanted to add multiple different sidechains, put those isotopically
labelled sidechains in a file and use this reaction
```
scaffold {
  smarts: "[cH]"
}
sidechain {
  smarts: "[1]"
  join { a1: 0 a2: 0}
}
```
This time the invocation might look like
```
consecutive_reactions -P file.rxn -N 2 file.smi sidechains.smi > new.smi
```
Make sure to make your query specific enough. The above is very generic. When used on
1000 random molecules, and with just 5 sidechains, it took 66 seconds and generated 556k products.

Use the -v option to see what options are being given to trxn. You can add more options
via `-trxn ... -trxn` of grab the default command generated here and configure to needs.

For example if you wished the second reaction to only add the second substituent to an atom in
the same ring system as the one onto which the first has been added, use the same reaction as
above for the first step. The second reacton might be:
```
scaffold {
  id: 0
  query {
    query {
      ring_system_specifier {
        base {
          environment: "[R]!@[1]"
          set_global_id: 1
        }
      }
      smarts: "[/IWgid1cH]"
    }
  }
}
sidechain {
  id: 1
  smarts: "[1]"
  join { a1: 0 a2: 0}
}
```
Provide this as the -P2 argument, which specifies a second reaction used for all
reactions after the first.
This time it takes 23 seconds to generate 232k molecules from the
starting 100 and 5 different sidechains. Beware of combinatorics.

If you wished to restrict the second substituent to the same ring, use `ring_specifier` rather
than the `ring_system_specifier` used above.

A word of explation on the reaction file.

Normally for such a complex query, it would be placed in a separate file and in the
proto specification we would use a `query_file: "fname"` construct. That also facilitates debugging
with tsubstructure. But it is also possible to inline the query as done here.

The ring system specifier looks for a ring system that contains a ring atom, bonded to an
atom with an isotope - something from the first reaction. All atoms in that ring system
are assigned global id 1.

The smarts then looks for an atom with global id 1 and which is an aromatic carbon
with an attached Hydrogen atom. This ensures that the new addition goes onto an atom
in the same ring system as where the first one was attached.

That invocation might look like
```
consecutive_reactions.sh -v -P add_sidechain.rxn -P2 add_to_same_ring.rxn \
        file.smi sidechain.smi
```

The following options are recognised;

-P <rxn>                name of a textproto reaction file
-P2 <rxn>               [optional] name of a second textproto reaction file.
                        Used for all reactions except the first.
-N <number>             number of consecutive invocations (default 2)
-trxn ... -trxn         passed to all invocations of trxn, eg: -trxn -v -W NONE -trxn 
-v                      verbose output
MSG

  $stderr << msg << "\n"
  exit(0)
end

def main
  cl = IWCmdline.new("-v-N=ipos-P=sfile-P2=sfile-trxn=close")
  if cl.unrecognised_options_encountered
    $stderr << "unrecognised_options_encountered\n";
    usage
  end

  verbose = cl.option_present('v')

  unless cl.option_present('P')
    $stderr << "must specify reaction via the -P option\n"
    usage
  end

  rxn = cl.value('P')

  rxn2 = cl.value('P2')
  $stderr << "RXN2 #{rxn2}\n"

  n = if cl.option_present('N')
        cl.value('N')
      else
        2
      end

  trxn = "trxn.sh -d -d -m each -z i -J wrscaf -J nomshmsg -S -"
  trxn1 = "#{trxn} -P #{rxn}"
  $stderr << "trxn1 #{trxn1}\n"
  if cl.option_present('P2')
    rxn2 = cl.value('P2')
    trxn2 = "#{trxn} -P #{rxn2}"
  else
    trxn2 = trxn1
  end
  $stderr << "trxn2 #{trxn2}\n"

  if cl.option_present('trxn')
    trxn << ' ' << cl.value('trxn')
  end

  scaffold = ARGV.shift
  sidechains = ARGV.join(' ')

  cmd = "#{trxn1} #{scaffold} #{sidechains}"
  (n - 1).times.each do 
    cmd << "| #{trxn2} -i smi - #{sidechains}"
  end
  $stderr << "Executing #{cmd}\n" if verbose
  system(cmd)
end

main
