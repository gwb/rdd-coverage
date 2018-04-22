module NYCPlots
    using GeoInterface
    import PyPlot; plt=PyPlot
    plt.rc("figure", dpi=300.0)
    # plt.rc("figure", figsize=(6,4))
    # plt.rc("figure", autolayout=true)
    plt.rc("savefig", dpi=300.0)
    plt.rc("text", usetex=true)
    plt.rc("font", family="serif")
    plt.rc("font", serif="Palatino")
    cbbPalette = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7"]
    using LaTeXStrings
    using LibGEOS
    using ..NYC: read_distr_reprpoints
    using Distributions: MultivariateNormal
    
    using Formatting
    
    using PyCall
    PyCall.unshift!(PyVector(pyimport("sys")["path"]), "")
    @PyCall.pyimport NYC_prepare_plots
    background_schdistrs = NYC_prepare_plots.background_schdistrs
    
    function plot_polygon(poly_coords, color, alpha)
        ax = plt.gca()
        poly_array = Array{Float64}(length(poly_coords), 2)
        for i in 1:length(poly_coords)
            poly_array[i, :] = poly_coords[i]
        end
        polygon = plt.plt[:Polygon](poly_array, true)
        p = ax[:add_patch](polygon)
        p[:set_color](color)
        p[:set_zorder](99)
        p[:set_alpha](alpha)
        p[:set_edgecolor]("none")
    end

    function plot_buffer(shape, border, color, dist, alpha)
        border_buffer = LibGEOS.buffer(border, dist)
        shape_buffer = LibGEOS.intersection(border_buffer, shape)            
        if typeof(shape_buffer) <: LibGEOS.Polygon
            poly_coords = GeoInterface.coordinates(shape_buffer)[1]
            plot_polygon(poly_coords, color, alpha)
        elseif typeof(shape_buffer) <: LibGEOS.MultiPolygon
            mpoly_coords = GeoInterface.coordinates(shape_buffer)
            for poly_coords in mpoly_coords
                plot_polygon(poly_coords[1], color, alpha)
            end
        else
            println(typeof(shape_buffer))
        end
    end
    
    function plot_schdistr_labels(;fontsize=12.0, labelcolor="black")
        plt.rc("text", usetex=false)
        schdistr_represent = read_distr_reprpoints()
        xmin,xmax = plt.xlim()
        ymin,ymax = plt.ylim()
        for schdistr in keys(schdistr_represent)
            center_x, center_y = schdistr_represent[schdistr]
            
            if !(xmin < center_x < xmax)
                continue
            elseif !(ymin < center_y < ymax)
                continue
            end
            
            plt.annotate(
                @sprintf("%d", schdistr),
                xy=(center_x, center_y),
                xycoords="data",
                xytext=(0,0),
                textcoords="offset points",
                horizontalalignment="center",
                verticalalignment="center",
                fontweight="black",
                color=labelcolor,
                alpha=1.0,
                fontsize = fontsize,
                zorder=3
            )
            plt.annotate(
                @sprintf("%d", schdistr),
                xy=(center_x, center_y),
                xycoords="data",
                xytext=(0,0),
                textcoords="offset points",
                horizontalalignment="center",
                verticalalignment="center",
                fontweight="black",
                color=labelcolor,
                alpha=0.5,
                fontsize = fontsize,
                zorder=20
            )
        end
        plt.rc("text", usetex=true)
    end

    function plot_all_pairs(τpost_pairs, sentinels, borders_dict, schdistr_shape_dict; scaleup=200.0)
        for distr_pair in keys(sentinels)
            distrA, distrB = distr_pair

            τpost = τpost_pairs[distr_pair]
            eff_size = abs(mean(τpost)) / std(τpost) * 2
            if eff_size < 2.0
                continue
            end

            if distrA ∈ (1,2) || distrB ∈ (1,2)
                # don't do Manhattan, it's weird
                continue
            end

            if distrA >= distrB
                continue
            end
            X◫ = hcat(sentinels[distr_pair].coordinates...)
            border = borders_dict[(parse(distrA), parse(distrB))]

            border_coords = GeoInterface.coordinates(border)
            if typeof(border) <: LibGEOS.MultiLineString
            else
                border_coords = [border_coords]
            end
            ax = plt.gca()

            distrA_shape = schdistr_shape_dict[parse(distrA)]
            distrB_shape = schdistr_shape_dict[parse(distrB)]

            if mean(τpost) < 0
                cheaper_shape = distrA_shape
                pricier_shape = distrB_shape
            else
                cheaper_shape = distrB_shape
                pricier_shape = distrA_shape
            end
            plot_buffer(pricier_shape, border, cbbPalette[2], (abs(mean(τpost)) + std(τpost))*scaleup, 0.6)
            plot_buffer(pricier_shape, border, cbbPalette[1], abs(mean(τpost))*scaleup, 0.6)
        end
    end

    function plot_streets(;zoom=11, brighten=30, rgb=(200, 200, 200), zorder=-9)
        ax = plt.gca()
        ax[:set_autoscale_on](false)
        ax[:set_autoscalex_on](false)
        ax[:set_autoscaley_on](false)
        EPSG=2263 # projection
        topleft =  (40.75, -74.01)
        botright = (40.58, -73.60)
        py"""
        import matplotlib.pyplot as plt

        import geopandas as gpd
        import tilemapbase
        import pyproj
        import PIL.Image as _Image

        def bounds(plotter, epsg):
            epsg_proj = pyproj.Proj(init='epsg:%d' % epsg, preserve_units=True)
            epsg_3857 = pyproj.Proj(init='epsg:%d' % tilemapbase.mapping._WEB_MERCATOR)
            self = plotter
            scale = 2 ** self.zoom
            x0, y0 = self.xtilemin / scale, self.ytilemin / scale # in tile space
            x1, y1 = (self.xtilemax + 1) / scale, (self.ytilemax + 1) / scale
            x0_ll, y0_ll = tilemapbase.to_lonlat(x0, y0)
            x1_ll, y1_ll = tilemapbase.to_lonlat(x1, y1)
            y0_ll, y1_ll = y1_ll, y0_ll # for some reason that I don't understand
            print(y0_ll)
            print(y1_ll)
            x0_epsg, y0_epsg = epsg_proj(x0_ll, y0_ll)
            x1_epsg, y1_epsg = epsg_proj(x1_ll, y1_ll)
            return x0_epsg, x1_epsg, y0_epsg, y1_epsg
            
        def brighten(img, factor):
            # split the image into individual bands
            source = img.split()
            R, G, B, A = 0, 1, 2, 3

            # process each band separately
            out_R = source[R].point(lambda x: min(x+factor, 256))
            out_G = source[G].point(lambda x: min(x+factor, 256))
            out_B = source[B].point(lambda x: min(x+factor, 256))

            # build a new multiband image
            brightened = _Image.merge(img.mode, (out_R, out_G, out_B, source[A]))
            return brightened

        def makecolor(img, rgb):
            # split the image into individual bands
            source = img.split()
            R, G, B, A = 0, 1, 2, 3

            # process each band separately
            out_R = source[R].point(lambda x: rgb[0])
            out_G = source[G].point(lambda x: rgb[1])
            out_B = source[B].point(lambda x: rgb[2])

            # build a new multiband image
            new_im = _Image.merge(img.mode, (out_R, out_G, out_B, source[A]))
            return new_im

        tilemapbase.start_logging()
        tilemapbase.init("/Users/imolk/tmp/tilemapbase.cache", create=True)

        extent = tilemapbase.Extent.from_lonlat($topleft[1], $botright[1],
                          $botright[0], $topleft[0])
        # plotter_ln = tilemapbase.Plotter(extent,tilemapbase.tiles.Carto_Dark_No_Labels, zoom=$zoom)
        plotter_ln = tilemapbase.Plotter(extent,tilemapbase.tiles.Stamen_Terrain_Lines, zoom=$zoom)
        plotter = plotter_ln
        # tile = brighten(plotter.as_one_image(), $brighten)
        tile = makecolor(plotter.as_one_image(), $rgb)
        $ax.imshow(tile, interpolation="lanczos", extent=bounds(plotter, $EPSG), origin="upper", zorder=$zorder)
        """
    end
    function plot_sales(filtered, missing_sqft; edgecolor="black", labelsize=12.0, labelcolor="black", colorbar=true)
        fig=plt.figure()
        background_schdistrs(plt.gca(),
                    color="#AAAAAA", 
                    edgecolor=edgecolor,
                    alpha=1.0,
                    linestyle="-",
                    zorder=-10
                    )
        background_schdistrs(plt.gca(),
                    color="none", 
                    edgecolor=edgecolor,
                    alpha=1.0,
                    linestyle="-",
                    linewidth=0.5,
                    zorder=11)
        _xlim = (0.98e6, 1.06e6)
        _ylim = (1.5e5, 2.2e5)
        fig[:set_size_inches](7.0,5.0)
        plt.plot(
            missing_sqft[:XCoord], 
            missing_sqft[:YCoord], 
            color="white",
            marker="x",
            markersize=1.0,
            linestyle="",
            zorder=8,
            alpha=0.7,
            linewidth=0.1,
            markeredgewidth=0.3,
        )
        plt.scatter(filtered[:XCoord], 
            filtered[:YCoord], 
            c=filtered[:logSalePricePerSQFT],
            cmap="jet",
            edgecolors="None",
            zorder=10,
            marker="o",
        #     vmin=4.0,
        #     vmax=8.0,
            alpha=1.0,
            s=1,
            )
        plt.xlim(_xlim)
        plt.ylim(_ylim)
        # plt.xlim(np.percentile(sales_geocoded_ggl["XCoord"].dropna(), [1,99]))
        # plt.ylim(np.percentile(sales_geocoded_ggl["YCoord"].dropna(), [1,99]))
        # plt.axis("off")
        # plt.title("log(Property Price per SQFT) in New York")
        ax = plt.gca()
        ax[:set](aspect="equal", adjustable="box")
        ax[:get_xaxis]()[:set_ticks]([])
        ax[:get_yaxis]()[:set_ticks]([])
        ax[:set_facecolor]("#B0DAEE")
        if colorbar
            cbar = plt.colorbar(label="Square foot prices (\\\$ per sqft)")
            ax = cbar[:ax]
            tick_template = [1, 2, 3, 4, 5, 6, 7, 8, 9]
            labeled = [10, 100, 1000, 10^4]
            ax[:set_yticks](log.(labeled))
            ax[:set_yticks](log.([tick_template; tick_template.*10; tick_template.*100; tick_template.*1000]), 
                                   minor=true)
            ax[:set_yticklabels]([format(x, commas=true) for x in labeled])
        end
        
        ax[:set_autoscale_on](false)
        ax[:set_autoscalex_on](false)
        ax[:set_autoscaley_on](false)
        plot_schdistr_labels(fontsize=labelsize, labelcolor=labelcolor)
    end
    function plot_cliffface(μ, Σ, color; label=L"posterior of $\tau(x)$ on residuals")
        plt.plot(μ, color=color, ".", label=label)
        msize=10
        plt.plot(1,μ[1], color=cbbPalette[5], 
            markersize=msize, marker="o")
        plt.plot(length(μ),μ[end], color=cbbPalette[4], 
            markersize=msize, marker="o")
        posterior_sd = sqrt.(diag(Σ))
        plt.fill_between(1:length(μ), μ.-2*posterior_sd, μ.+2*posterior_sd, 
            color=color, 
            alpha=0.3,
            linewidth=0,
            zorder=-2,
            )
        posterior_distr = MultivariateNormal(μ, Σ)
        for _ in 1:10
            plt.plot(1:length(μ), rand(posterior_distr), "-", color="white", alpha=0.4, linewidth=1.5, zorder=-1)
        end
        plt.xlim(1,length(μ))
        plt.ylim(minimum(μ.-2*posterior_sd), maximum(μ.+2*posterior_sd))
        
        # second y-axis for percentage price increase
        ax = plt.gca()
        y1, y2=ax[:get_ylim]()
        x1, x2=ax[:get_xlim]()
        ax2=ax[:twinx]()
        ax2[:set_ylim](y1, y2)
        yticks_frac = collect(0.5: 0.1 : 2.0)
        yticks_log = log.(yticks_frac)
        in_window = y1 .< yticks_log .< y2
        ax2[:set_yticks](yticks_log[in_window])
        ylabels = [@sprintf("%.0f\\%%", y*100) for y in yticks_frac[in_window]]
        ax2[:set_yticklabels](ylabels)
        
        plt.sca(ax)
    end
end
