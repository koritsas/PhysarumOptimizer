import org.apache.commons.math3.analysis.UnivariateFunction;

/**
 * Created by ilias on 16/9/2016.
 */
public class DiameterChoosingFunction implements UnivariateFunction {

    public double value(double D) {
        boolean negative=false;
        if (D<0){
            negative=true;
        }

        D=Math.abs(D);

        if (D<0.05){
            D=0;
        }else if(D<0.110){
            D=0.110;
        }else if(D<0.14){
            D=0.14;
        }else if(D<0.16){
            D=0.16;
        }else if(D<0.2){
            D=0.2;
        }else if(D<0.225){
            D=0.225;
        }else if(D<0.28){
            D=0.28;
        }else if(D<0.315){
            D=0.315;
        }else if(D<0.355){
            D=0.355;
        }else if(D<0.4){
            D=0.4;
        }else if(D<0.5){
            D=0.5;
        }else{
            D=0.6;
        }

        if(negative){
            D=-D;
        }

        return D;
    }
}
