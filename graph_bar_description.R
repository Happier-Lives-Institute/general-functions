# A bar graph 
# Neat bars across the different scale outcomes
# Can show the means
# Can show the counts

# Set colours for the graph if there are no grouping conditions
color_1 <- "#1C5FB8"
color_2 <- "#ff6900"

bar_description <- function(
        data, variable,
        grouping = "",
        xlabel = "",
        flip = F, meanToLeft = T,
        individual_labels = T,
        max_x = 10
){
    # Need the data as data.frame for certain manipulations
    data <- as.data.frame(data)

    # Calculate text distance based on max_x
    text_distance <- 0.5 * max_x/10

    # Prepare the graph (according to whether there's a grouping variable or not)
    if(grouping == ""){
        p <- data %>% ggplot(aes(x=get(variable))) +
            geom_bar(fill=color_1, color=color_1)
    } else {
        p <- data %>% ggplot(aes(
            x=get(variable),
            color = get(grouping),
            fill = get(grouping)
        )) +
            geom_bar(position = position_dodge2(preserve = "single")) +
            labs(fill = grouping, color = grouping)
    }

    p <- p + xlab(xlabel) +
        cowplot::theme_cowplot()

    # If this is a numeric variable set the x axis and add the mean
    if(is.numeric(data[,variable])){

        # Depending on whether there is grouping get information
        if(grouping == ""){
            # Get max count
            max_count <- max((data %>% count(get(variable)))$n)

            # Get the mean
            my_mean <- mean(data[,variable], na.rm=T)
        } else {
            # Get max count (of the highest condition)
            counts <- data %>% count(get(grouping), get(variable))
            max_count <- max(counts$n)

            # Get the mean
            my_mean <- data %>%
                group_by(get(grouping)) %>%
                summarise(mean=mean(get(variable), na.rm=T), n = n()) %>%
                rowwise() %>%
                mutate(
                    text_position = ifelse(meanToLeft, mean - text_distance, mean + text_distance),
                    label = round_c(mean, 2),
                ) %>% ungroup()
            names(my_mean)[1] <- grouping
        }

        # Set up the scaling
        p <- p +
            scale_y_continuous(limits = c(0, max_count + 5)) +
            scale_x_continuous(breaks = 0:max_x,
                               limits = c(0, max_x),
                               # limits = c(-0.1, 10.1),
                               # To make sure the bars at 0 and 10 appear
                               oob = scales::rescale_none
            )


        # Depending on grouping, set the means
        if(grouping == ""){
            p <- p +
                geom_vline(
                    xintercept = my_mean,
                    linetype="dashed", color=color_2, linewidth=1
                ) +
                geom_text(
                    x = ifelse(meanToLeft, my_mean - text_distance, my_mean + text_distance),
                    y = max_count+3,
                    label=round_c(my_mean),
                    color = color_2
                )
        } else {
            p <- p +
                geom_vline(
                    data = my_mean,
                    aes(xintercept = mean, color = get(grouping)),
                    linetype="dashed", linewidth=1
                ) +
                geom_text(
                    data = my_mean,
                    aes(
                        x = text_position,
                        label=label,
                        color = get(grouping)
                    ),
                    y = max_count+3
                )
        }
    }

    # Should there be individual column labels?
    if(individual_labels) {
        # Set the count for individual columns according to whether the plot has
        # been flipped or not.
        if(flip){
            p <- p + coord_flip() + geom_text(
                stat='count',
                aes(label=..count..),
                hjust= +1.05,
                color="white"
            )
        } else {
            p <- p + geom_text(stat='count', aes(label=..count..), vjust=-0.25, color="black")
        }
    }

    return(p)
}
