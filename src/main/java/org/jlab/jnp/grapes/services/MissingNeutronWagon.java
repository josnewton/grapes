package org.jlab.jnp.grapes.services;

import org.jlab.jnp.hipo4.data.Bank;
import org.jlab.jnp.hipo4.data.Event;
import org.jlab.jnp.hipo4.data.SchemaFactory;
import java.util.ArrayList;
import java.util.HashMap;

/**
 * 
 * Skim for studying high-energy pion rejection for improving electron identification 
 *
 * Until standard filters provide access to momentum and ECAL energy.
 *
 * @author jnewton
 * @author baltzell
 */
public class MissingNeutronWagon extends Wagon {

    static final double BEAM_ENERGY = 10.6f;
    static final double PROTON_MASS = 0.938f;
    static final double ELE_MOM_LOW = 0.5f;
    static final double ELE_MOM_HIGH = 4.5f;
    static final double PION_MOM_LOW = 5.0f;
    static final double NEUTRON_MASS_LOW = 0.7f;
    static final double NEUTRON_MASS_HIGH = 1.3f;
    
    public MissingNeutronWagon() {
        super("MissingNeutronWagon","jnewton","0.2");
    }

    @Override
    public boolean init(String jsonString) {
        System.out.println("MissingNeutronWagon READY.");
        return true;
    }

    private double getMomentum(final int ipart, Bank particles) {
        final double px = particles.getFloat("px",ipart);
        final double py = particles.getFloat("py",ipart);
        final double pz = particles.getFloat("pz",ipart);
        return (double)Math.sqrt(px*px + py*py + pz*pz);
    }

    double getMissingMass(final int ipart1, final int ipart2, Bank particles) {
	 
        final double px1 = particles.getFloat("px",ipart1);
        final double py1 = particles.getFloat("py",ipart1);
        final double pz1 = particles.getFloat("pz",ipart1);
        final double e1 = Math.sqrt(px1*px1 + py1*py1 + pz1*pz1);

        final double px2 = particles.getFloat("px",ipart2);
        final double py2 = particles.getFloat("py",ipart2);
        final double pz2 = particles.getFloat("pz",ipart2);
        final double e2 = Math.sqrt(px2*px2 + py2*py2 + pz2*pz2);

        double missing_energy = BEAM_ENERGY + PROTON_MASS - (e1 + e2);
        double missing_px = px1 + px2;
        double missing_py = py1 + py2;
        double missing_pz = BEAM_ENERGY - (pz1 + pz2);
        double missing_p = Math.sqrt(missing_px*missing_px + missing_py*missing_py + missing_pz*missing_pz);
        double missing_mass = Math.sqrt(missing_energy*missing_energy - missing_p*missing_p);
        
        return missing_mass;
    }

    @Override
    public boolean processDataEvent(Event event, SchemaFactory factory) {

        Bank particles = new Bank(factory.getSchema("REC::Particle"));
        event.read(particles);
        if (particles.getRows()<1) return false;
        
        ArrayList<Integer> eleCandi = new ArrayList<>();
        ArrayList<Integer> posCandi = new ArrayList<>();
        ArrayList<Integer> negCandi = new ArrayList<>();

        for (int ipart=0; ipart<particles.getRows(); ipart++) {
           
            final int pid=particles.getInt("pid",ipart);
            final int status=particles.getInt("status",ipart);
            final boolean isFD = (int)(Math.abs(status)/1000) == 2;

            // count electron/positive/negative candidates based on EB pid:
            if (isFD) {
                if      (pid==11) eleCandi.add(ipart);
                else if (pid==211 || pid==-11) posCandi.add(ipart);

                // FIXME: this first row assumption is broken, trigger particle can be in any row.
                else if (pid==-211 || (pid==11 && ipart>0)) negCandi.add(ipart);
            }
        }

        // abort asap:
        if (eleCandi.isEmpty() && posCandi.isEmpty() && negCandi.isEmpty()) return false;

        ///////////////////////////////////////////////////////////////////
        ///////////////////////////////////////////////////////////////////
  
        // Positive Pion:
        if (eleCandi.size()>0 && posCandi.size()>0) return true;
	
        //Negative Pion:
        if (eleCandi.size()>0 && negCandi.size()>0) return true;

        // FIXME:  why do anything else after this point?
        // Everything below is already covered by the previous returns!
        
        // electrons, negatives, and positives within momentum cuts:
        ArrayList<Integer> eleHiCandi = new ArrayList<>();
        ArrayList<Integer> negHiCandi = new ArrayList<>();
  
        for (int ii : eleCandi) {
            if (this.getMomentum(ii,particles) > ELE_MOM_LOW && this.getMomentum(ii,particles) < ELE_MOM_HIGH) {
                eleHiCandi.add(ii);
            }
        }
        for (int ii : negCandi) {
            if (this.getMomentum(ii,particles) > PION_MOM_LOW) negHiCandi.add(ii);
        }

        // e-, negative
        if (eleHiCandi.size()>0 && negHiCandi.size()>0) return true;
	
        ArrayList<Integer> posHiNeutron = new ArrayList<>();
	     
        for(int ii : eleCandi) {
            for(int jj : posCandi) {
                final double mm = this.getMissingMass(ii,jj,particles);
                if (mm>NEUTRON_MASS_LOW && mm<NEUTRON_MASS_HIGH) {
                    posHiNeutron.add(jj);
                }
            }
        }

        //electron, positive track, and correct missing mass 
        if(eleHiCandi.size()>0 && posHiNeutron.size()>0) return true;
       
        return false;
    }

}
