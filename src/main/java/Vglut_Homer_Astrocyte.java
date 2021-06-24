/*
 * Compute distances from Vglut-Homer synapses to astrocyte border
 * Author Philippe Mailly
 */

import fiji.util.gui.GenericDialogPlus;
import ij.*;
import ij.gui.Plot;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.PlugIn;
import ij.plugin.filter.RankFilters;
import ij.process.ImageConverter;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.awt.Polygon;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.ThreadLocalRandom;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Objects3DPopulation;
import mcib3d.geom.Point3D;
import mcib3d.geom.Voxel3D;
import mcib3d.image3d.ImageFloat;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import mcib3d.image3d.distanceMap3d.EDT;
import mcib3d.image3d.processing.MaximaFinder;
import mcib3d.utils.ArrayUtil;



public class Vglut_Homer_Astrocyte implements PlugIn {

    private boolean canceled = false;
    private String imageDir = "";
    private String outDirResults = "";
    private final Calibration cal = new Calibration();
    
// min distance synapse
    private double minDistSynap = 0.4;
// do median filter
    private boolean median = false;
// max distance to find dots to astrocyte
    private final double maxDist = 1.5;
// min intensity of vglutSted dot in vglut confocal channel
    private  double vglutConfDotsIntRef = 10.00;
// min intensity of homerSted dot in homer confocal channel
    private  double homerConfDotsIntRef = 6.00;
// min distance to defie a synapse
    private final double synDist = 0.3;   
    private double maxTolerance = 15;   
    private final ArrayList<Double> vglutHomerDist = new ArrayList();
    private BufferedWriter outPutResults;
    public double VglutTol, HomerTol;
    private float rxy = 2.5f;
    private float rz = 2.5f;
    private final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));


    private boolean dialog() {
        canceled = true;
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 180, 0);
        gd.addImage(icon);
        gd.addDirectoryField("Choose images directory", imageDir);
        gd.addCheckbox(" Median filter", median);
        gd.addNumericField("Vglut confocal channel reference intensity : ", vglutConfDotsIntRef, 2);
        gd.addNumericField("Homer confocal channel reference intensity : ", homerConfDotsIntRef, 2);
        gd.addNumericField("Minimum distance between 2 synapses : ", minDistSynap, 3);
        gd.showDialog();
        if (gd.wasCanceled())
            canceled = false;
        imageDir = gd.getNextString()+ File.separator;
        median = gd.getNextBoolean();
        vglutConfDotsIntRef = gd.getNextNumber();
        homerConfDotsIntRef = gd.getNextNumber();
        minDistSynap = gd.getNextNumber();
        return(canceled);
    }
    
    
 /**
  * Find dots sted terminaisons
  * take only dots sted if intensity in confocal channel > dotsIntRef 
  */
    private void findDotsIntPop (ImagePlus imgConf, ArrayList<Point3D> points, double intRef) {
        if (points.size() > 0) {
            for (int n = 0; n < points.size(); n++) {
                int x = (int)points.get(n).x;
                int y = (int)points.get(n).y;
                double confInt = imgConf.getProcessor().getPixelValue(x, y);
                if (confInt < intRef) {
                    points.remove(n);
                    n--;
                }
            }
        }
    }
    
  
    /**
     * return objects population in an binary image
     * @param img
     * @return pop objects population
     */

    private static  Objects3DPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        Objects3DPopulation pop = new Objects3DPopulation(labels);
        return pop;
    }   
     
    // Flush and close images
    private void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    }

    
     /*Median filter 
     * 
     * @param img
     * @param size
     */ 
    public static void median_filter(ImagePlus img, double size) {
        RankFilters median = new RankFilters();
        for (int s = 1; s <= img.getNSlices(); s++) {
            img.setZ(s);
            median.rank(img.getProcessor(), size, RankFilters.MEDIAN);
            img.updateAndDraw();
        }
    }
    
    
    /** 
     * Find sted dots
     * @param img channel
     * @return ArrayList<Point3D> stedDots
     */
    public ArrayList<Point3D> findStedDots(ImagePlus img, String maxTol, int imageNum) {
        if (median)
            median_filter(img, 0.5);
        ArrayList<Point3D> pt3Ds = new ArrayList();
        
        if (img.getNSlices() == 1) {
            Roi roiPt;
            if (imageNum == 1) {
                img.show();
                IJ.run(img, "Find Maxima...", "");
                double tol = IJ.getNumber(maxTol+" tolerance = ", 0);
                img.hide();
                if (maxTol.equals("vglut"))
                    VglutTol = tol;
                else
                    HomerTol = tol; 
            }
            else {
                if (maxTol.equals("vglut"))
                    maxTolerance = VglutTol;
                else
                    maxTolerance = HomerTol;
                IJ.run(img, "Find Maxima...", "prominence="+maxTolerance+" exclude output=[Point Selection]");    
            }
            roiPt = img.getRoi();
            Polygon poly = roiPt.getPolygon();
            for (int n = 0; n < poly.npoints; n++) {
                double x = (double)poly.xpoints[n];
                double y = (double)poly.ypoints[n];
                pt3Ds.add(new Point3D(x, y, 1));   
            }
            FileSaver imgMax = new FileSaver(img);
            imgMax.saveAsTiff(outDirResults+img.getTitle()+"_Max.tif");    
        }
        else {
            if (imageNum == 1) {
                img.show();
                img.setSlice(img.getNSlices()/2);
                IJ.run(img, "Find Maxima...", "");
                double tol = IJ.getNumber(maxTol+" Tolerance = ", 0);
                img.hide();
                if (maxTol.equals("vglut"))
                    VglutTol = tol;
                else
                    HomerTol = tol; 
                maxTolerance = tol;
                img.deleteRoi();
            }
            else {
                if (maxTol.equals("vglut"))
                    maxTolerance = VglutTol;
                else
                    maxTolerance = HomerTol;
            }
            ImagePlus imgColor = img.duplicate();
            imgColor.setCalibration(cal);
            IJ.run(imgColor, "16-bit", "");
            ImageConverter imgConv = new ImageConverter(imgColor);
            imgConv.convertToRGB();
            ImageHandler imgH = ImageHandler.wrap(img);
            imgH.setCalibration(cal);
            MaximaFinder test = new MaximaFinder(imgH, (float)maxTolerance);
            test.setRadii(rxy, rz);
            test.setVerbose(false);
            //test.getImagePeaks().show();
            // list
            ArrayList<Voxel3D> list = test.getListPeaks();  
            imgH.closeImagePlus();
            // draw point
            for (Voxel3D vox : list) {
                int x = vox.getRoundX();
                int y = vox.getRoundY();
                int z = vox.getRoundZ()+1;
                Point3D pt = new Point3D(x, y, z);
                pt3Ds.add(pt);
                imgColor.setSlice(z);
                imgColor.updateAndDraw();
                ImageProcessor ip = imgColor.getProcessor();
                ip.setColor(Color.yellow);
                ip.drawOval(x-3, y-3, 6, 6);
            }
            FileSaver imgMax = new FileSaver(imgColor);
            imgMax.saveAsTiff(outDirResults+img.getTitle()+"_Max.tif");
            flush_close(imgColor);
        }      
        return(pt3Ds);
    } 
    
    /**
     * Astrocyte Cell
     * 
     */
    private void detectAstro(ImagePlus img) {
        String stack = "";
        if (img.getNSlices() > 1) {
            stack = " stack";
            img.setSlice(img.getNSlices()/2);
        }
        IJ.run(img, "16-bit", "");
        IJ.run(img, "Difference of Gaussians", "  sigma1=0.30 sigma2=0.05 scaled"+stack);
        IJ.setAutoThreshold(img, "Li dark"); 
        Prefs.blackBackground = false;
        IJ.run(img, "Convert to Mask", "method=Li background=Dark");
    }
    
    
    /**
     * keep synapses at distance sup minDistSynap
     */
    
    private ArrayList<Point3D> removeSynapes(ArrayList<Point3D> synPop) {
        ArrayList<Point3D> synPt = new ArrayList<>();
        for (int n = 0; n < synPop.size(); n++) {
            boolean distOk = false;
            Point3D ptA = synPop.get(n);
            for (int i = 0; i < synPop.size(); i++) {
                if (n != i) {
                    Point3D ptB = synPop.get(i);
                    double dist = ptA.distance(ptB, cal.pixelWidth, cal.pixelDepth);
                    if (dist > minDistSynap) {
                        distOk = true;
                    }
                    else {
                        distOk = false;
                        i = synPop.size();
                    }
                }
            }
            if (distOk)
                synPt.add(ptA);
        }
        return(synPt);
    }
    
    
     /**
     * Find synapses 
     * if distance between vglut dots and homer dots <= synDist add to synapPop
     * else remove vglut and homer object
     * keep only objects that make synapses
     * @param vglutPop
     * @param homerPop
     * @return 
     */
    
    private ArrayList<Point3D> findSynapses (ArrayList<Point3D> vglutPoints, ArrayList<Point3D> homerPoints) {
        double dist = 0;
        int homerIndex = 0;
        int vglutIndex = 0;
        Point3D homerPt = new Point3D();
        Point3D vglutPt = new Point3D();
        ArrayList<Point3D> synapsesPop = new ArrayList();
        // find synapses from smaller pop
        IJ.showStatus("Finding synapses ...");
        if (homerPoints.size() < vglutPoints.size()) { 
            for (int i = 0; i < homerPoints.size(); i++) {
                double distMin = Double.MAX_VALUE;
                homerPt = homerPoints.get(i);
                // Find nearest vglut Point from homer Point
                for (int n = 0; n < vglutPoints.size(); n++) {
                    vglutPt = vglutPoints.get(n);
                    dist = homerPt.distance(vglutPt, cal.pixelWidth, cal.pixelDepth);
                    if (dist < distMin) {
                        distMin = dist;
                        homerIndex = i;
                        vglutIndex = n;
                    }
                       
                }
                // if dist <= synapse Dist compute center distance point
                if (distMin <= synDist) {
                    Point3D synapse = new Point3D((vglutPoints.get(vglutIndex).x + homerPoints.get(homerIndex).x)/2, (vglutPoints.get(vglutIndex).y + homerPoints.get(homerIndex).y)/2,
                            (vglutPoints.get(vglutIndex).z + homerPoints.get(homerIndex).z)/2);
                    synapsesPop.add(synapse);
                    vglutHomerDist.add(distMin);
                    System.out.println("Distance = "+dist);
                }
            }
        }
        else {
            for (int i = 0; i < vglutPoints.size(); i++) {
                double distMin = Double.MAX_VALUE;
                vglutPt = vglutPoints.get(i);
                // Find nearest vglut Point from homer Point
                for (int n = 0; n < homerPoints.size(); n++) {
                    homerPt = homerPoints.get(n);
                    dist = vglutPt.distance(homerPt, cal.pixelWidth, cal.pixelDepth);
                    if (dist < distMin) {
                       distMin = dist;
                       homerIndex = n;
                       vglutIndex = i;
                    }
                }
                // if dist <= synapse Dist compute center distance point
                if (distMin <= synDist) {
                    Point3D synapse = new Point3D((vglutPoints.get(vglutIndex).x + homerPoints.get(homerIndex).x)/2, (vglutPoints.get(vglutIndex).y + 
                            homerPoints.get(homerIndex).y)/2, (vglutPoints.get(vglutIndex).z + homerPoints.get(homerIndex).z)/2);
                    synapsesPop.add(synapse);
                    vglutHomerDist.add(distMin);
                    //System.out.println("Distance = "+dist);
                }
            }
        }
        System.out.println(synapsesPop.size()+" synpases found"); 
        ArrayList<Point3D> synDistOk = removeSynapes(synapsesPop);
        System.out.println(synDistOk.size()+" synpases found after distance filter"); 

    return(synDistOk);
    }
        
    
    
    /** Find min distance from synapses dots to astro border using EDT
     * 
     *  Compute Volumes of synapses dots
     * @param synapse dots population
     * @return results ArrayList of measures synapse dot volume, distance ....
     */
    private ArrayList<Double> distanceSynToAstroBorder (ImagePlus imgAstro, Objects3DPopulation astroPop, ArrayList<Point3D> synapsePoints) {
        double dist;
        ArrayList<Double> synapseDist = new ArrayList();
        // EDT inverse astro info 
        ImageHandler img = ImageInt.wrap(imgAstro).createSameDimensions();
        astroPop.draw(img, 128);
        ImageFloat edtAstro = EDT.run(img, 120, (float)cal.pixelWidth, (float)cal.pixelDepth, true, 0);
        
        for (int i = 0; i < synapsePoints.size(); i++) { 
            IJ.showStatus("Finding distances ...");
            Point3D synapsePoint = synapsePoints.get(i);
            // distance to astro
            int x = (int)(synapsePoint.x);
            int y = (int)(synapsePoint.y);
            int z = (int)(synapsePoint.z);
            edtAstro.getImagePlus().setSlice(z);

            dist = edtAstro.getImagePlus().getProcessor().getPixelValue(x, y);
            //System.out.println((i+1)+"x= "+x+" y= "+y+" z= "+z+" dist = "+dist);
            // if dist > distMax remove synapse
            if (dist > maxDist)  {
                synapsePoints.remove(synapsePoint);
                i--;
            }
            else {
                synapseDist.add(dist);
            }
        }
        
//        FileSaver edtSave = new FileSaver(edtAstro.getImagePlus());
//        edtSave.saveAsTiff(outDirResults + img.getTitle() + "_EDT" + ".tif");
        edtAstro.closeImagePlus();
        img.closeImagePlus();
        System.out.println(synapsePoints.size()+" synapses arround "+maxDist);
        return(synapseDist);
    }
    
    
    /**
     * Find synapses closest distances 
     */
    
    private ArrayUtil closestDistances(ArrayList<Point3D> synPt) {
        ArrayUtil synDist = new ArrayUtil(synPt.size());
        for (int n = 0; n < synPt.size(); n++) {
            double closeDist = Double.MAX_VALUE;
            Point3D ptA = synPt.get(n);
            for (int i = 0; i < synPt.size(); i++) {
                if (n != i) {
                    Point3D ptB = synPt.get(i);
                    double dist = ptA.distance(ptB, cal.pixelWidth, cal.pixelDepth); 
                    if (dist < closeDist) {
                        closeDist = dist;
                    }
                }
            }
            synDist.addValue(n, closeDist);
        }
        return(synDist);
    }
    
    /**
     * Get distances statistics from Synapses
     */
    
    private ArrayUtil distStats(ArrayList<Point3D> synapsesPt, String imgName) {
        ArrayUtil dist = closestDistances(synapsesPt);
        double mean = dist.getMean();
        double sd = dist.getStdDev();
        DecimalFormat dec = new DecimalFormat("#0.0000");
        ArrayUtil[] histo = dist.getHistogram(100);
        Plot dPlot = new Plot("Distances plot", "Closest distances", "");
        dPlot.setLimits(0, histo[0].getMaximum(), 0, histo[1].getMaximum());
        dPlot.setColor(Color.black, Color.black);
        dPlot.add("bars", histo[0].getArray(), histo[1].getArray());
        dPlot.addLabel(0.4, 0.2, "Mean = "+ dec.format(mean) + " +- " + dec.format(sd));
        dPlot.draw();
        ImagePlus imgPlot = dPlot.getImagePlus();
        FileSaver plotSave = new FileSaver(imgPlot);
        plotSave.saveAsTiff(outDirResults + imgName + "_ClosestDist_Plot" + ".tif");
        flush_close(imgPlot);  
        return(dist);
    }
    
    
    private void savePopImage(ImagePlus img, ArrayList<Point3D> synapsePoint) {
        String imageName = img.getTitle();
        ImagePlus imgColor = img.duplicate();
        IJ.run(imgColor, "16-bit", "");
        ImageConverter imgConv = new ImageConverter(imgColor);
        imgConv.convertToRGB();
        // draw point
        for (int i = 0; i < synapsePoint.size(); i++) {
            imgColor.setSlice(synapsePoint.get(i).getRoundZ());
            ImageProcessor ip = imgColor.getProcessor();
            ip.setColor(Color.yellow);
            int x = synapsePoint.get(i).getRoundX();
            int y = synapsePoint.get(i).getRoundY();
            ip.drawOval(x-2, y-2, 4, 4);
        }
        FileSaver ObjectsFile = new FileSaver(imgColor);
        ObjectsFile.saveAsTiff(outDirResults +File.separator+ imageName + "_Mask.tif");
        imgColor.close();
    }
    
        
   // tag obj number with random color
    private void tagsObject(ImagePlus img, ArrayList<Point3D> synapses) {  
        Font tagFont = new Font(Font.MONOSPACED, Font.PLAIN, 8);
        int col;
        for (int i = 0; i < synapses.size(); i++) {
            Point3D obj = synapses.get(i);
            col = ThreadLocalRandom.current().nextInt(1, 255);
            int x = obj.getRoundX();
            int y = obj.getRoundY();
            int z = obj.getRoundZ()+1;
            img.setZ(z);
            // draw point and label
            ImageProcessor ip = img.getProcessor();
            ip.setColor(col);
//            ip.drawOval(x, y, 2, 2);
            ip.drawDot(x, y);
            ip.setFont(tagFont);
            ip.drawString(Integer.toString((i+1)), (x-4), (y));
            img.updateAndDraw();
        }    
    }
    
    /**
     * Save objects Population in image
     * @param astroPop
     * @param synapsePop
     * @param imageName
     * @param img 
     */
    private void saveImgObjects(Objects3DPopulation astroPop, ArrayList<Point3D> synapsePoints, String imageName, ImagePlus img) {
        //create image objects population
        ImageHandler imgObj = ImageInt.wrap(img).createSameDimensions();
        imgObj.set332RGBLut();
        imgObj.setCalibration(cal);
        
        //Astro population gray
        if (astroPop.getNbObjects() > 0)
            astroPop.draw(imgObj, 255);

        //Synapses population        
        if (synapsePoints.size() > 0)
            tagsObject(imgObj.getImagePlus(), synapsePoints);
        
        //save image
        FileSaver ObjectsFile = new FileSaver(imgObj.getImagePlus());
        if (imgObj.getImagePlus().getNSlices() == 1)
            ObjectsFile.saveAsTiff(outDirResults +File.separator+ imageName + "_Objects.tif");
        else
            ObjectsFile.saveAsTiffStack(outDirResults +File.separator+ imageName + "_Objects.tif");

        /*
        * close images
        */
        imgObj.closeImagePlus();
    }
    
    /** Write results
     * 
     */
    private void writeResults(String image, Objects3DPopulation astroPop, ArrayList<Point3D> vglutPop, ArrayList<Point3D> homerPop,
                 ArrayList<Double> synapseDist, ArrayUtil distStats) throws IOException {
        FileWriter fileResults = null;
        
        try {
            fileResults = new FileWriter(outDirResults + image +"_results.xls", false);
        } catch (IOException ex) {
            Logger.getLogger(Vglut_Homer_Astrocyte.class.getName()).log(Level.SEVERE, null, ex);
        }
        outPutResults = new BufferedWriter(fileResults);
        outPutResults.write("Astrocyte Vol\tSynapse distance Vglut-Homer\tSynapse Dist to Astro border\tSynapse Closest dist\n");
        outPutResults.flush();
        int[] maxObj = {astroPop.getNbObjects(), vglutPop.size(), homerPop.size()};
        ArrayUtil max = new ArrayUtil(maxObj);
        int  nbObj = (int)max.getMaximum();
        for (int n = 0; n < nbObj; n++) {
            if (n < astroPop.getNbObjects())
                outPutResults.write(astroPop.getObject(n).getVolumeUnit()+"\t");
            else
                 outPutResults.write("-\t");
            if (n < vglutHomerDist.size())
                outPutResults.write(vglutHomerDist.get(n)+"\t");
            else 
                outPutResults.write("-\t");
            if (n < synapseDist.size())
                outPutResults.write(synapseDist.get(n)+"\t");
            else 
                outPutResults.write("-\t");
            if (n < distStats.getSize())
                outPutResults.write(distStats.getValue(n)+"\n");
            else 
                outPutResults.write("-\n");
        }
        outPutResults.flush();
        outPutResults.close();
    }
    
    /**
     * 
     * @param arg
     */
    @Override
    public void run(String arg) {
        if (!dialog()) {
            IJ.showMessage(" Pluging canceled");
            return;
        }

        File inDir = new File(imageDir);
        String[] imageFile = inDir.list();
        if (imageFile == null) {
            return;
        }
        // create output folder
        outDirResults = inDir + File.separator+ "Out"+ File.separator;
        File outDir = new File(outDirResults);
        if (!Files.exists(Paths.get(outDirResults))) {
            outDir.mkdir();
        }
        
        Arrays.sort(imageFile);
        int imageNum = 0;
        for (String f : imageFile) {
            // open image files
            if (f.endsWith(".ics")) {
                imageNum++;
                try {
                    // for all image files
                    // get file name without extension
                    String rootName = f.replace(".ics", "");
                    // Find spatial calibration
                    // create OME-XML metadata store of the latest schema version
                    ServiceFactory factory;
                    factory = new ServiceFactory();
                    OMEXMLService service = factory.getInstance(OMEXMLService.class);
                    IMetadata meta = service.createOMEXMLMetadata();
                    ImageProcessorReader reader = new ImageProcessorReader();
                    reader.setMetadataStore(meta);
                    // read calibration metadata in first file
                    String imageName = imageDir+f;
                    reader.setId(imageName);
                    double sx = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
                    double sy = meta.getPixelsPhysicalSizeY(0).value().doubleValue();
                    double sz= meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
                    cal.pixelWidth = sx;
                    cal.pixelHeight = sy;
                    cal.pixelDepth = sz;
                    cal.setUnit("microns");
                    System.out.println("Image "+rootName);
                    System.out.println(" x ="+cal.pixelWidth+" z ="+cal.pixelDepth);
                    
                    /* read channels
                    *  C0 = Astro 488
                    *  C1 = Vglut 647 Confocal
                    *  C2 = Vglut 647 STED
                    *  C3 = Homer 594 Confocal
                    *  C4 = Homer 594 STED
                    */
                    // Open channels
                    // Astro channel
                    int series = 0;
                    reader.setSeries(series);
                    ImporterOptions options = new ImporterOptions();
                    options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                    options.setId(imageName);
                    options.setSplitChannels(true);
                    
                    options.setCBegin(series, 0);
                    options.setCEnd(series, 0);
                    ImagePlus imgC0 = BF.openImagePlus(options)[0];
                    imgC0.setTitle("Astro");
                    IJ.showStatus("Finding astrocyte cell ...");
                    detectAstro(imgC0);
                    // Find Astro pop
                    Objects3DPopulation astroPop = getPopFromImage(imgC0);
                    System.out.println("Astro ="+astroPop.getNbObjects());

                    // Vglut
                    options.setCBegin(series, 1);
                    options.setCEnd(series, 2);
                    ImagePlus imgC1 = BF.openImagePlus(options)[0];
                    imgC1.setTitle("VglutConf");
                    ImagePlus imgC2 = BF.openImagePlus(options)[1];
                    imgC2.setTitle("VglutSted");
                    IJ.showStatus("Finding vglut spots ...");
                    ArrayList<Point3D> vglutStedPoints = findStedDots(imgC2, "vglut", imageNum);
                    System.out.println("Vglut ="+vglutStedPoints.size());
                    
                    // take only vglutSted dots if intensity in vglutConf image > dotsIntRef
                    findDotsIntPop(imgC1, vglutStedPoints, vglutConfDotsIntRef);
                    System.out.println("Vglut in confocal ="+vglutStedPoints.size());
                    // save VglutStedPop on VglutConf image
                    savePopImage(imgC1,vglutStedPoints);

                    flush_close(imgC1);
                    flush_close(imgC2);

                    // Homer
                    options.setCBegin(series, 3);
                    options.setCEnd(series, 4);
                    ImagePlus imgC3 = BF.openImagePlus(options)[0];  
                    imgC3.setTitle("HomerConf");
                    ImagePlus imgC4 = BF.openImagePlus(options)[1];
                    imgC4.setTitle("HomerSted");
                    IJ.showStatus("Finding Homer spots ...");
                    ArrayList<Point3D> homerStedPoints = findStedDots(imgC4, "homer", imageNum);
                    System.out.println("Homer ="+homerStedPoints.size());
                    
                    
                    // take only HomerSted dots if intensity in homerConf image > dotsIntRef
                    findDotsIntPop(imgC3, homerStedPoints, homerConfDotsIntRef);
                    System.out.println("Homer in confocal ="+homerStedPoints.size());
                    // save HomerStedPop on HomerConf image
                    savePopImage(imgC3,homerStedPoints);
                    flush_close(imgC3);
                    flush_close(imgC4);
                    
                    // Find synapses
                    // synapse = distance homerSted dots vglutSted dots <= synapDist
                    IJ.showStatus("Finding synapses ....");
                    ArrayList<Point3D> synapsesPoint= findSynapses(vglutStedPoints, homerStedPoints);
                    
                    // Compute distance synapses to border astro
                    IJ.showStatus("Computing distances ....");
                    // with Distance map
                    ArrayList<Double> dist = distanceSynToAstroBorder(imgC0, astroPop, synapsesPoint);
                    

                    // Compute statistical distribution for synapses pop
                    ArrayUtil distStats = new ArrayUtil(synapsesPoint.size());
                    if (synapsesPoint.size() > 1)
                        distStats = distStats(synapsesPoint, rootName);
                    
                    try {
                        // Writing results
                        writeResults(rootName, astroPop, vglutStedPoints, homerStedPoints, dist, distStats);
                    } catch (IOException ex) {
                        Logger.getLogger(Vglut_Homer_Astrocyte.class.getName()).log(Level.SEVERE, null, ex);
                    }
                    // Save image objects
                    saveImgObjects(astroPop, synapsesPoint, rootName, imgC0);
                    flush_close(imgC0);
                    vglutHomerDist.clear();
                    } catch (DependencyException | ServiceException | FormatException | IOException ex) {
                        Logger.getLogger(Vglut_Homer_Astrocyte.class.getName()).log(Level.SEVERE, null, ex);
                    }
            }
        }
        IJ.showStatus("Process done");
    }    
}
